#' Initialize tile scores for aggregation
#' 
#' Higher scores favor merging.
#' 
#' Methods for different modes of aggregation:
#'  1. `dscore` has already been manually calculated and stored in `aggs$edges$dscore`.
#'     Will set the `score`, `score_size`, and `compactness` attributes of `aggs$meta_data` to 0.
#'  2. `dscore` is the product of three factors:
#'     * `w`: A 2-cluster GMM is used to determine the mean `mu` and standard deviation `sig` of the distance
#'       between tiles that have similar gene expression (Euclidean distance `d` in PC space).
#'       Then we define `d_mu = mu + sig` and `d_sig = alpha * sig`, and calculate `w` as
#'       `w = 0.5 - 1 / (1 + exp(-(d - d_mu) / d_sig))`. Ranges from -0.5 to 0.5. If adjacent tiles
#'       are very dissimilar (`d >> d_mu`), then `d` is large, and `w` is close to `-0.5`. If adjacent
#'       tiles are very similar (`d < d_mu`), then `d` is small, and `w` is positive. 
#'     * `score_size`: `(1 - npts_from/max_npts) * (1 - npts_to/max_npts)`. Ranges from 0 to 1.
#'       If merging the two tiles would have a total number of points ≥`max_npts`, then `score_size` is `-Inf`,
#'       which prevents merging.
#'     * `dC`: `.5 * (C_merge - C_from - C_to + 1)`. Ranges from 0 to 1.
#'  3. `dscore` is the product of same three factors as mode 2, but `score_size` has already been precomputed.
#'     Additionally, `dscore` is set to -1 if both adjacent tiles have at least `min_npts` cells.
#' 
#' @param aggs A tile data structure.
#' @param agg_mode Method to use to calculate aggregation scores (1, 2, or 3).
#' @param ... Additional parameters for different modes:
#'  1. No additional parameters.
#'  2. Requires `alpha` and `max_npts`. For `alpha`, 0.2 = conservative merging, 2 = liberal merging.
#'  3. Requires `alpha` and `min_npts`. For `alpha`, 0.2 = conservative merging, 2 = liberal merging.
#' 
#' @returns `aggs` with scores for aggregation stored in attributes:
#'   \item{edges}{Additional attributes are calculated:
#'     * `dscore`: Overall score for merging two tiles. Product of `w`, `score_size`, and `dC`.
#'     * `w`: Gene expression similarity score.
#'     * `score_size`: Penalizes tiles with many points.
#'     * `perimeter_merge`: Perimeter of merged tile.
#'   }
#'   \item{d_mu,d_sig}{Parameters used to calculate `w`.}
#'   \item{pcs_merged}{Average PCs for merged tile.}
#' 
#' @export
init_scores = function(aggs, agg_mode, ...) {
    ## (1) Gene expression Scores
    # aggs$edges$ll = as.numeric(get_edges_ll(
    #     aggs$pcs, 
    #     aggs$meta_data$npts, 
    #     aggs$edges$from-1,
    #     aggs$edges$to-1, 
    #     fname
    # ))
    # aggs$edges$entropy = as.numeric(get_edges_entropy(
    #     aggs$pcs, 
    #     aggs$meta_data$npts, 
    #     aggs$edges$from-1,
    #     aggs$edges$to-1
    # ))

    ## OLD: based on cosine similarity between clusters 
    # aggs$pcs_norm = sweep(aggs$pcs, 1, sqrt(rowSums(aggs$pcs^2)), '/')
    # X = sweep(aggs$pcs, 1, sqrt(rowSums(aggs$pcs^2)), '/')
    # aggs$edges$score_cosine = rowSums(aggs$pcs_norm[aggs$edges$from, ] * aggs$pcs_norm[aggs$edges$to, ]) ## -1 (far apart) to 1 (identical)
    # aggs$edges$score_cosine = 1 - .5 * (rowSums((aggs$pcs_norm[aggs$edges$from, ] - aggs$pcs_norm[aggs$edges$to, ])^2)) ## -1 (far apart) to 1 (identical)

    if (agg_mode == 1) {
        if (!'dscore' %in% colnames(aggs$edges)) {
            stop('for agg_mode=1, dscore must already be initialized')
        }
        ## set up dummy variables 
        aggs$meta_data$score = 0
        aggs$meta_data$score_size = 0
        aggs$meta_data$compactness = 0
        
    } else if (agg_mode == 2) {
        alpha = list(...)[['alpha']]
        max_npts = list(...)[['max_npts']]
        stopifnot(!is.null(alpha))
        stopifnot(!is.null(max_npts))
        
        
        # GMM to decide which distances are good or bad 
        d = sqrt(rowSums((aggs$pcs[aggs$edges$from, ] - aggs$pcs[aggs$edges$to, ])^2))
        mres = mclust::Mclust(d, G=2)
        ord = order(mres$parameters$mean)
        stopifnot(mres$parameters$mean[ord[1]] < mres$parameters$mean[ord[2]])
        d_sig = sqrt(mres$parameters$variance$sigmasq)[ord[1]]
        d_mu = mres$parameters$mean[ord[1]] + d_sig
        d_sig = alpha * d_sig


        ## w is designed so that all d < d_mu will be kept and all d > d_mu will be cut. 
        ## w := probability that edge is a good internal edge (between -0.5 [large d] and 0.5 [small d])
        aggs$edges$w = 1 / (1 + exp(-(d - d_mu) / d_sig)) 
        aggs$edges$w = 0.5 - aggs$edges$w ## high distance = low probability 

        ## For the merged scores, compute the same thing
        ## Use the same cached function to turn distances into weights 

        ## Weighted sum of PCs
        a0 = aggs$meta_data$npts[aggs$edges$from]
        a1 = aggs$meta_data$npts[aggs$edges$to]
        atot = a0 + a1
        aggs$pcs_merged = sweep(aggs$pcs[aggs$edges$from, ], 1, a0, '*') + sweep(aggs$pcs[aggs$edges$to, ], 1, a1, '*')
        aggs$pcs_merged = sweep(aggs$pcs_merged, 1, atot, '/')
        aggs$edges$perimeter_merge = aggs$meta_data$perimeter[aggs$edges$from] + aggs$meta_data$perimeter[aggs$edges$to] - 2 * aggs$edges$edge_length
        # aggs$edges$nedges_internal_merge = 1


        f_size_from = (-1/max_npts) * aggs$meta_data$npts[aggs$edges$from] + 1
        f_size_to = (-1/max_npts) * aggs$meta_data$npts[aggs$edges$to] + 1
        aggs$edges$score_size = f_size_from * f_size_to


        ## Change in compactness
        C_from = (a0 / atot) * (4 * pi * aggs$meta_data$area[aggs$edges$from]) / (aggs$meta_data$perimeter[aggs$edges$from]^2) 
        C_to = (a1 / atot) * (4 * pi * aggs$meta_data$area[aggs$edges$to]) / (aggs$meta_data$perimeter[aggs$edges$to]^2)
        C_merge = (4 * pi * aggs$edges$area) / (aggs$edges$perimeter_merge ^ 2)
        dC = .5 * (C_merge - C_from - C_to + 1) ## ranges from 0 to 1

        aggs$d_mu = d_mu
        aggs$d_sig = d_sig
        aggs$edges$dscore = aggs$edges$w * aggs$edges$score_size * dC
        # set dscore to -Inf if merge agg exceeds size thresholds
        # (-1 is lower than the lowest possible dscore, and allows the edge to be preserved in output)
        aggs$edges$dscore[aggs$edges$npts >= max_npts] = -1
    } else if (agg_mode == 3) {
        alpha = list(...)[['alpha']]
        min_npts = list(...)[['min_npts']]
        stopifnot(!is.null(alpha))
        stopifnot(!is.null(min_npts))
        
        
        ## First, re-initialize all dscore values 
        a0 = aggs$meta_data$npts[aggs$edges$from]
        a1 = aggs$meta_data$npts[aggs$edges$to]
        atot = a0 + a1

        # GMM to decide which distances are good or bad 
        d = sqrt(rowSums((aggs$pcs[aggs$edges$from, ] - aggs$pcs[aggs$edges$to, ])^2))
        mres = mclust::Mclust(d, G=2)
        ord = order(mres$parameters$mean)
        stopifnot(mres$parameters$mean[ord[1]] < mres$parameters$mean[ord[2]])
        d_sig = sqrt(mres$parameters$variance$sigmasq)[ord[1]]
        d_sig = alpha * d_sig
        d_mu = mres$parameters$mean[ord[1]] + d_sig
        aggs$edges$w = 1 / (1 + exp(-(d - d_mu) / d_sig)) 
        aggs$edges$w = 1 - aggs$edges$w ## high distance = low probability 

        ## Score_size should already be correctly set 
        ## ... do nothing ... 

        ## Change in compactness
        C_from = (a0 / atot) * (4 * pi * aggs$meta_data$area[aggs$edges$from]) / (aggs$meta_data$perimeter[aggs$edges$from]^2) 
        C_to = (a1 / atot) * (4 * pi * aggs$meta_data$area[aggs$edges$from]) / (aggs$meta_data$perimeter[aggs$edges$from]^2)
        C_merge = (4 * pi * aggs$edges$area) / (aggs$edges$perimeter_merge ^ 2)
        dC = .5 * (C_merge - C_from - C_to + 1) ## ranges from 0 to 1

        aggs$edges$dscore = aggs$edges$w * aggs$edges$score_size * dC

        aggs$d_mu = d_mu
        aggs$d_sig = d_sig
        ## Delete edges between aggs of good size 
        # (-1 is lower than the lowest possible dscore, and allows the edge to be preserved in output)
        aggs$edges$dscore[aggs$meta_data$npts[aggs$edges$from] >= min_npts & aggs$meta_data$npts[aggs$edges$to] >= min_npts] = -1
    } else {
        stop(glue('agg mode {agg_mode} not implemented yet'))
    }
    
    
#     ## (2) Size Scores
#     aggs$meta_data$score_size = dnorm(aggs$meta_data$npts, mean = size_mu, sd = size_mu, log = TRUE) - size_factor
#     aggs$edges$score_size = dnorm(aggs$edges$npts, mean = size_mu, sd = size_mu, log = TRUE) - size_factor


#     ## (3) Compactness Scores
#     aggs$meta_data$compactness = (4 * pi * aggs$meta_data$area) / (aggs$meta_data$perimeter * aggs$meta_data$perimeter) 
#     perim = aggs$meta_data$perimeter[aggs$edges$from] + aggs$meta_data$perimeter[aggs$edges$to] - 2 * aggs$edges$edge_length
#     # perim = aggs$meta_data$perimeter[aggs$edges$from] + aggs$meta_data$perimeter[aggs$edges$to] - 2 * st_length(st_sfc(aggs$edges$shape))
#     aggs$edges$compactness = (4 * pi * aggs$edges$area) / (perim * perim) 


#     ## (4) Combine parts of scores 
#     ## No exp needed, do everything in log space
#     ## First for aggs 
#     aggs$meta_data$score = weight_size * aggs$meta_data$score_size + weight_compactness * aggs$meta_data$compactness + weight_entropy * aggs$meta_data$entropy
#     ## Then for edges, in 3 parts: 
#     ##   1. score of merged ij
#     ##   2. subtract individual scores of i and j 
#     ##   3. add joint composition score 
#     aggs$edges$score = weight_size * aggs$edges$score_size + weight_compactness * aggs$edges$compactness + weight_entropy * aggs$edges$entropy

#     ## Compute the difference in LL if we merge 
#     ## Weight scores by number of cells 
#     # aggs$edges$
#     score_from = aggs$meta_data$score[aggs$edges$from]
#     # aggs$edges$
#     score_to = aggs$meta_data$score[aggs$edges$to]
#     w = cbind(
#         aggs$meta_data$npts[aggs$edges$from],
#         aggs$meta_data$npts[aggs$edges$to]
#     ) %>% prop.table(1)
#     # aggs$edges$dscore = aggs$edges$score - rowSums(aggs$edges[, .(score_from, score_to)] * w)        
#     aggs$edges$dscore = aggs$edges$score - rowSums(cbind(score_from, score_to) * w) 
    
    # ## UNWEIGHTED
    # aggs$edges$dscore = aggs$edges$score - aggs$meta_data$score[aggs$edges$from] - aggs$meta_data$score[aggs$edges$to]

    # ## WEIGHTED
    # alpha = prop.table(cbind(aggs$meta_data$npts[aggs$edges$from], aggs$meta_data$npts[aggs$edges$to]), 1)    
    # aggs$edges$dscore = aggs$edges$score - alpha[, 1] * aggs$meta_data$score[aggs$edges$from] - alpha[, 2] * aggs$meta_data$score[aggs$edges$to]
    # aggs$edges$dscore = aggs$edges$dscore + weight_cosine * aggs$edges$score_cosine 
    
    return(aggs)   
}


# update_agg = function(aggs, e_merge) {
#     ## Update PCs
#     alpha = prop.table(aggs$meta_data$npts[c(e_merge$from, e_merge$to)])    
#     aggs$pcs[e_merge$from, ] = alpha[1] * aggs$pcs[e_merge$from, ] + alpha[2] * aggs$pcs[e_merge$to, ]
#     aggs$pcs_norm[e_merge$from, ] = aggs$pcs[e_merge$from, ] / sqrt(sum(aggs$pcs[e_merge$from, ]^2))
    
#     ## Update size score 
#     aggs$meta_data$score_size[e_merge$from] = dnorm(aggs$meta_data$npts[e_merge$from], mean = size_mu, sd = size_mu, log = TRUE) - size_factor

#     ## Update compactness score 
#     perim = aggs$meta_data$perimeter[e_merge$from] + aggs$meta_data$perimeter[e_merge$to] - 2 * e_merge$edge_length 
#     aggs$meta_data$compactness[e_merge$from] = (4 * pi * e_merge$area) / (perim * perim)    

#     ## Summarize agg score
#     aggs$meta_data$score[e_merge$from] = weight_size * aggs$meta_data$score_size[e_merge$from] + 
#                                          weight_compactness * aggs$meta_data$compactness[e_merge$from]
    
#     return(aggs)
# }

# update_scores = function(aggs, e_update) {
#     ## (1) Composition Scores 
#     # X = sweep(aggs$pcs, 1, sqrt(rowSums(aggs$pcs^2)), '/')

#     aggs$edges$score_cosine[e_update] = 1 - .5 * (Matrix::rowSums((aggs$pcs_norm[aggs$edges$from[e_update], ] - aggs$pcs_norm[aggs$edges$to[e_update], ])^2)) ## -1 (far apart) to 1 (identical)

#     ## (2) Size Scores
#     aggs$edges$score_size[e_update] = dnorm(aggs$edges$npts[e_update], mean = size_mu, sd = size_mu, log = TRUE) - size_factor

#     ## (3) Compactness Scores
#     perim = aggs$meta_data$perimeter[aggs$edges$from[e_update]] + aggs$meta_data$perimeter[aggs$edges$to[e_update]] - 2 * aggs$edges$edge_length[e_update] #st_length(st_sfc(aggs$edges$shape[e_update]))
#     aggs$edges$compactness[e_update] = (4 * pi * aggs$edges$area[e_update]) / (perim * perim) 
    
#     ## (4) Combine parts of scores 
#     ## No exp needed, do everything in log space
#     ## Then for edges, in 3 parts: 
#     ##   1. score of merged ij
#     ##   2. subtract individual scores of i and j 
#     ##   3. add joint composition score 
#     aggs$edges$score[e_update] = weight_size * aggs$edges$score_size[e_update] + weight_compactness * aggs$edges$compactness[e_update] ## already affinity, so no negative sign 

#     # ## UNWEIGHTED
#     # aggs$edges$dscore[e_update] = aggs$edges$score[e_update] - aggs$meta_data$score[aggs$edges$from[e_update]] - aggs$meta_data$score[aggs$edges$to[e_update]]
    
#     ## WEIGHTED 
#     alpha = prop.table(cbind(aggs$meta_data$npts[aggs$edges$from[e_update]], aggs$meta_data$npts[aggs$edges$to[e_update]]), 1)    
#     aggs$edges$dscore[e_update] = aggs$edges$score[e_update] - alpha[, 1] * aggs$meta_data$score[aggs$edges$from[e_update]] - alpha[, 2] * aggs$meta_data$score[aggs$edges$to[e_update]]
    
#     aggs$edges$dscore[e_update] = aggs$edges$dscore[e_update] + weight_cosine * aggs$edges$score_cosine[e_update] 
    
#     return(aggs)    
# }




#' @export
pcs_to_rgb = function(V) {
    rgb = uwot::umap(as.matrix(V), min_dist=.05, spread=.2, n_components = 3L, fast_sgd = TRUE)
    rgb = sweep(rgb, 2, apply(rgb, 2, min), '-') 
    rgb = sweep(rgb, 2, apply(rgb, 2, max), '/')     
    rgb = round(255 * rgb)
    rgb = as.matrix(rgb)
    rgb = round(t(rgb2hsv(rgb[, 1], rgb[, 2], rgb[, 3])) * 255) ## hsv
    rgb = rgb2hex(rgb)
    return(rgb)
    
}

#' @export
color_aggs = function(aggs, seed=1) {
    set.seed(seed)
    ## First, color discrete components with different colors
    ##        no two adjacent aggs should have the same color
    # g = igraph::from_data_frame()$fun(dplyr::select(aggs$edges, from, to), directed = FALSE, )
    g = igraph::from_data_frame()$fun(
        rbind(
            dplyr::select(aggs$edges, from, to), ## add self edges to not remove disconnected aggs
            cbind(from = 1:nrow(aggs$meta_data), to = 1:nrow(aggs$meta_data))
        ),
        directed = FALSE
    )
    y = igraph::greedy_vertex_coloring(g)
    y = y[order(as.integer(names(y)))]
    while (length(unique(y)) < min(20, length(y))) {
        id_replace = names(which.max(table(y)))
        i = which(y == id_replace)
        y[i] = max(y) + sample(2, length(i), TRUE)
    
    }
    y = as.integer(factor(y))
    aggs$meta_data$color = tableau_color_pal('Classic 20')(length(unique(y)))[y]
    
    ## Next, color by continuous scale from gene expression 
    aggs$logcpx = normalizeData(aggs$counts, median(Matrix::colSums(aggs$counts)), 'log')
    aggs$pca = weighted_pca(aggs$logcpx, rep(1, nrow(aggs$meta_data)), npc = 10, do_corr = TRUE)    
    aggs$meta_data$rgb = pcs_to_rgb(aggs$pca$embeddings)
    aggs$logcpx = NULL
    aggs$pca = NULL

    ## Also, if you already have PC embeddings (from agglom step), get those colors too 
    if ('pcs' %in% names(aggs)) {
        aggs$meta_data$rgb2 = pcs_to_rgb(aggs$pcs)
    }
    
    return(aggs)
}

#' Merges tiles using single-linkage agglomerative clustering
#' 
#' Note that this function mutates many of the input values to keep
#' track of updated values for tiles and their borders after
#' each successive merge. Every merge of adjacent tiles has an
#' associated score (higher means merging is more favorable).
#' Merging is conducted greedily, one step of a time, updating
#' the score associated with pair of tiles at each step.
#' 
#' @returns The `aggs` data structure with updated values for
#'   \item{meta_data}{A data table with attributes for the new merged tiles:
#'     * `X,Y`: Centroid of each tile. (NOT updated)
#'     * `npts`: Number of points in each tile. (Updated)
#'     * `shape`: A `sfc` list-column with the geometries for each tile. (NOT updated)
#'     * `area,perimeter`: Area and perimeter of each tile. (Updated)}
#'   \item{edges}{A data table with attributes for edges between adjacent merged tiles,
#'     where edges with infinite `dscore` are removed:
#'     * `from,to`: Tile IDs for the two tiles bordering this edge. (Updated for new merged tiles)
#'     * `x0,y0,x1,y1`: Centroid coordinates for the two tiles bordering this edge. (NOT updated)
#'     * `area,npts`: Sum of areas and numbers of points in the two tiles bordering this edge. (Updated)
#'     * `edge_length`: Total length of the border between the `from` and `to` tiles. (Updated)
#'     * `dscore`: Overall score for merging two tiles. Product of `w`, `score_size`, and `dC`. (Updated)
#'     * `w`: Gene expression similarity score. (Updated)
#'     * `score_size`: Penalizes tiles with many points. (Updated)
#'     * `perimeter_merge`: Perimeter of merged tile. (Updated)}
#'   \item{pcs}{A `num_tiles` x `npcs` matrix with the average embedding value over all cells
#'     in each tile. (Updated)}
#'   \item{pcs_merged}{A `num_edges` x `npcs` matrix with average PCs for merged tiles. (Updated)}
#'   \item{d_mu,d_sig}{Parameters used to calculate `w`.}
#'   \item{aggmap}{A length `orig_num_tiles` vector mapping each original tile ID to the new
#'     tile IDs after merging.}
#' @export
merge_aggs = function(
    aggs, agg_mode, 
    d_mu=NULL, d_sig=NULL, iter_max=NULL,
    dscore_thresh=0,
    min_npts=0,
    max_npts=Inf,
    min_area=0, 
    max_area=Inf
) {
    ## without this, npts not updated by reference in C++ code
    ## b/c R int vector is cast as arma::vec and allocated different memory
    aggs$meta_data$npts = as.numeric(aggs$meta_data$npts) 
    aggs$edges$npts = as.numeric(aggs$edges$npts) 
    # aggs$meta_data$nedges = as.numeric(aggs$meta_data$nedges) 
    # aggs$edges$nedges = as.numeric(aggs$edges$nedges) 
    
    
    if (is.null(iter_max)) {
        iter_max = nrow(aggs$meta_data) - 1
    } else {
        iter_max = min(nrow(aggs$meta_data) - 1, iter_max)        
    }
    
    ## memory is a list of lists 
    ##   disjoint sets 
    ##   that results from iterative aggregation    
    memory = merge_aggs_cpp(
        V_pcs = aggs$pcs, 
        V_area = aggs$meta_data$area,
        V_perimeter = aggs$meta_data$perimeter,
        V_npts = aggs$meta_data$npts,
        # V_score = aggs$meta_data$score,
        # V_score_size = aggs$meta_data$score_size, 
        # V_compactness = aggs$meta_data$compactness, 
        # V_f_min = aggs$meta_data$f_min, 
        # V_nedges = aggs$meta_data$nedges, 
        # V_nedges_internal = aggs$meta_data$nedges_internal, 
        # V_f_dual_total = aggs$meta_data$f_dual_total, 

        # V_wl_ext = aggs$meta_data$wl_ext, 
        # V_w_ext = aggs$meta_data$w_ext, 
        # V_w_int = aggs$meta_data$w_int, 
        
        
        E_from = aggs$edges$from - 1,
        E_to = aggs$edges$to - 1,
        E_npts = aggs$edges$npts, 
        E_area = aggs$edges$area, 
        E_edge_length = aggs$edges$edge_length,
        # E_f_dual_total = aggs$edges$f_dual_total, 
        # E_nedges = aggs$edges$nedges, 
        # E_f_min_merge = aggs$edges$f_min_merge, 
        # E_nedges_merge = aggs$edges$nedges_merge, 
        # E_nedges_internal_merge = aggs$edges$nedges_internal_merge, 
        # E_f_dual_total_merge = aggs$edges$f_dual_total_merge, 
        
        E_pcs_merge = aggs$pcs_merge, 
        E_w = aggs$edges$w, 
        E_perimeter_merge = aggs$edges$perimeter_merge, 
        # E_wl_ext_merge = aggs$edges$wl_ext_merge, 
        # E_w_ext_merge = aggs$edges$w_ext_merge, 
        # E_w_int_merge = aggs$edges$w_int_merge, 
        E_score_size = aggs$edges$score_size,
        
        # E_score_merge = aggs$edges$score_merge,
        E_dscore = aggs$edges$dscore,
        
        d_mu=aggs$d_mu, ## for PCA based scoring 
        d_sig=aggs$d_sig, ## for PCA based scoring 
        # weight_size,
        # weight_compactness,
        # size_mu,
        # size_factor, 
        iter_max, 
        agg_mode,
        dscore_thresh,
        # min_area, 
        # max_area, 
        min_npts, 
        max_npts
    )
    aggs_keep = which(purrr::map_int(memory, length) > 0)
    memory = memory[aggs_keep]
    memory = purrr::map(memory, `+`, 1)
    
    o = order(Reduce(c, memory))
    aggmap = rep(seq_len(length(memory)), purrr::map_int(memory, length))[o]                               
    aggs$meta_data = aggs$meta_data[aggs_keep, ]
    aggs$meta_data$id = seq_len(nrow(aggs$meta_data))

    ## save this for mapping to points
    aggs$aggmap = aggmap
    
    e_keep = which(!is.infinite(aggs$edges$dscore))
    aggs$edges = aggs$edges[e_keep, ]
    aggs$pcs_merged = aggs$pcs_merged[e_keep, ]
    aggs$edges$from = aggmap[aggs$edges$from]
    aggs$edges$to = aggmap[aggs$edges$to]
    aggs$pcs = aggs$pcs[aggs_keep, ]
    
    return(aggs)

}

#' Update tile IDs in the DMT data structure after merging
#' 
#' @param dmt DMT data structure
#' @param agg Tile data structure after merging, with `aggmap` attribute
#'   mapping original tile IDs to new IDs for the merged tiles.
#' 
#' @returns DMT data structure with updated values for these attributes:
#'   * `pts$agg_id`: Unique ID for the tile that each point belongs to.
#'   * `edges$agg_from`: Unique ID for the tile that `edges$from_pt` belongs to.
#'   * `edges$agg_to`: Unique ID for the tile that `edges$to_pt` belongs to.
#' @export
update_dmt_aggid = function(dmt, aggs) {    
    dmt$pts$agg_id = aggs$aggmap[dmt$pts$agg_id]
    dmt$edges$agg_from = dmt$pts$agg_id[dmt$edges$from_pt]
    dmt$edges$agg_to = dmt$pts$agg_id[dmt$edges$to_pt]
    return(dmt)
}

#' Update shapes and counts matrix for tiles after merging
#' 
#' @param dmt DMT data structure
#' @param agg Tile data structure after merging
#' 
#' @returns `aggs` tile data structure with updated values for:
#'   \item{meta_data$shape}{Retraces the polygons for tiles using the new point assignments.}
#'   \item{meta_data$X,meta_data$Y}{Recomputes the centroid of each tile from the polygon shape.}
#'   \item{edges}{Ensures that `aggs$edges` always has `from` < `to`.}
#'   \item{counts}{A `num_genes` x `num_tiles` gene-by-tile matrix of aggregated transcript counts.}
#' @export
update_agg_shapes = function(dmt, aggs) {
    aggs$meta_data$shape = trace_polygons(dmt, aggs)
    xy = sf::st_coordinates(sf::st_centroid(aggs$meta_data$shape))
    aggs$meta_data$X = xy[, 'X']
    aggs$meta_data$Y = xy[, 'Y']
    aggs$edges[from > to, `:=`(from=to, to=from, x0=x1, y0=y1, x1=x0, y1=y0)]
    
    d = Matrix::sparse.model.matrix(~0 + factor(agg_id), dmt$pts)
    aggs$counts = dmt$counts %*% d
    colnames(aggs$counts) = gsub('^.*?(\\d+)$', '\\1', colnames(aggs$counts))
    
    return(aggs)
}
