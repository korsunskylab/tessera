#' Generic function that runs the Tessera algorithm on single-cell spatial data
#'
#' GetTiles is a generic function that runs the main Tessera algorithm.
#' If working with a Seurat object, please refer to the documentation of
#' the appropriate generic API: [GetTiles.Seurat()]. If users work with other
#' forms of the input, they can pass them directly to Tessera using the
#' [GetTiles.default()] API. The function arguments listed here are common in all
#' GetTiles interfaces.
#'
#' @family GetTiles
#' @rdname GetTiles
#' @inheritDotParams GetTiles.default -X -Y -counts -embeddings -loadings -meta_data
#'
#' @return If used with a Seurat object, it will return a pair of Seurat objects:
#' 1) the input single-cell object updated with tile assignments for each cell, and
#' 2) a Seurat object where each item represents an individual Tessera tile.
#'
#' For standalone operation, it returns Lists with the output of Tessera segmentation
#' (see [GetTiles.default()]).
#'
#' @export
GetTiles = function(...) {
    UseMethod("GetTiles")
}

#' Applies Tessera on a Seurat object
#'
#' @family GetTiles
#' @rdname GetTiles.Seurat
#'
#' @param obj Seurat object with spatial coordinates (and optionally, pre-computed
#'   single cell embeddings) stored as dimensional reductions.
#' @param spatial Name of dimensional reduction where the cells' x/y coordinates are stored.
#' @param embeddings Name of dimensional reduction where pre-computed single-cell embeddings are stored
#'   (a `num_cells` x `num_dim` matrix of cell embeddings across all latent dimensions).
#'   If missing, cell embeddings are calculated using PCA. If provided, the `npcs` parameter is ignored.
#' @param assay Seurat assay to pull data for when using the cell counts. Defaults to the DefaultAssay.
#' @param group.by Name of column in `obj@meta.data` to use for grouping cells into separate samples.
#' @param raw_results Whether to return the raw results from [GetTiles.default()].
#' @param tile.id.name Name of variable to store the tile IDs in the cell-level Seurat object.
#' @param reduction.name Name of dimensional reduction to store the aggregated tile-level embeddings
#'   in the tile-level Seurat object.
#' @param graph.name Name of graph to store tile adjacency matrix in the tile-level Seurat object.
#' @param add.isolated.cells Whether to add back isolated single cells that were pruned out.
#'   Only applies when `embeddings` are provided. Defaults to TRUE.
#' @inheritDotParams GetTiles.default -X -Y -counts -embeddings -loadings -meta_data
#'
#' @returns A List containing a pair of Seurat objects:
#' 1) `obj`: the input single-cell object whose meta.data has been updated with tile assignments for each cell
#' 2) `tile_obj`: a Seurat object where each item represents an individual Tessera tile
#'
#' @export
GetTiles.Seurat = function(
    obj,
    spatial,
    embeddings = NULL,
    dims.use = NULL,
    assay = NULL,
    group.by = NULL,
    raw_results = FALSE,
    tile.id.name = 'tile_id',
    reduction.name = 'pca',
    graph.name = 'tile_adj',
    add.isolated.cells = TRUE,
    ...
) {
    rlang::check_installed("Seurat", reason = "to use `GetTiles.Seurat()`")

    if ((!is.null(embeddings)) && (embeddings != 'harmony')) {
        warning('It is recommended to use harmonized cell embeddings as input for Tessera')
    }

    if (is.null(assay)) {
        assay = Seurat::DefaultAssay(obj)
    }

    if (!is.null(embeddings)) {
        emb <- Seurat::Embeddings(obj, reduction = embeddings)
        load <- Seurat::Loadings(obj, reduction = embeddings)
        
        if (is.null(dims.use)) {
            dims.use = seq_len(ncol(emb))
        }
        emb = emb[,dims.use,drop=FALSE]
        load = load[,dims.use,drop=FALSE]
    } else {
        emb <- NULL
        load <- NULL
    }

    if (is.null(group.by)) {
        warning('No value for group.by provided. Analyzing as a single sample.')
        if (is.null(obj@meta.data)) {
            group.by = 'group'
            obj@meta.data = data.frame(group = factor(rep(1, nrow(obj@meta.data))))
        } else {
            group.by = tail(make.unique(c(colnames(obj@meta.data), 'group')), n = 1) # avoid overwriting existing columns
            obj@meta.data[[group.by]] = factor(rep(1, nrow(obj@meta.data)))
        }
    } else {
        if (!(group.by %in% colnames(obj@meta.data))) {
            stop('group.by must be a column of meta.data')
        }
    }

    # Run Tessera
    res = GetTiles(
        X = Seurat::Embeddings(obj, spatial)[,1],
        Y = Seurat::Embeddings(obj, spatial)[,2],
        counts = obj[[assay]]$counts,
        embeddings = emb,
        loadings = load,
        meta_data = obj@meta.data,
        group.by = group.by,
        ...
    )

    if (raw_results) {return(res)}

    stopifnot(all(colnames(res$aggs$counts) == res$aggs$meta_data$id))

    # Collect tile metadata
    meta.data = data.frame(dplyr::select(res$aggs$meta_data, -shape))
    tile.shapes = res$aggs$meta_data$shape
    tile.counts = res$aggs$counts
    tile.embeddings = res$aggs$pcs
    row.names(meta.data) = meta.data$id
    stopifnot(all(colnames(res$aggs$counts) == res$aggs$meta_data$id))
    stopifnot(all(colnames(res$aggs$counts) == row.names(meta.data)))

    # Add tile IDs to input Seurat object
    if (tile.id.name %in% colnames(obj@meta.data)) {
        tile.id.name = tail(make.unique(c(colnames(obj@meta.data), tile.id.name)), n = 1) # avoid overwriting existing columns
        warning(paste0('To avoid overwriting existing meta.data variables, tile.id.name is set to ', tile.id.name))
    }
    obj@meta.data[[tile.id.name]] = as.character(NA)
    obj@meta.data[[tile.id.name]][res$dmt$pts$ORIG_ID] = res$dmt$pts$agg_id

    # (Optional) Add back isolated single cells that were pruned out
    if (add.isolated.cells & !is.null(embeddings)) {
        isolated_cells = which(is.na(obj@meta.data[[tile.id.name]]))

        if (length(isolated_cells) > 0) {
            warning(paste0('Adding back ', length(isolated_cells), ' isolated cells that were pruned out'))

            new_tiles_metadata = data.frame(
                id = seq_len(length(isolated_cells)) + ncol(tile.counts),
                X = Seurat::Embeddings(obj, spatial)[,1][isolated_cells],
                Y = Seurat::Embeddings(obj, spatial)[,2][isolated_cells],
                npts = 1,
                area = as.numeric(NA),
                perimeter = as.numeric(NA)
            )
            new_tiles_metadata[[group.by]] = obj@meta.data[[group.by]][isolated_cells]
            new_tiles_metadata$id = tail(
                make.unique(c(res$aggs$meta_data$id, as.character(new_tiles_metadata$id))),
                n = length(new_tiles_metadata$id)
            )
            rownames(new_tiles_metadata) = new_tiles_metadata$id

            obj@meta.data[[tile.id.name]][isolated_cells] = new_tiles_metadata$id

            new_counts = obj[[assay]]$counts[,isolated_cells,drop=FALSE]
            colnames(new_counts) = new_tiles_metadata$id
            stopifnot(all(rownames(new_counts) == rownames(tile.counts)))

            new_tiles_embeddings = emb[isolated_cells,]
            new_tiles_embeddings = do.call(cbind, replicate(  # concatenate embeddings to match tiles
                ncol(tile.embeddings) / ncol(new_tiles_embeddings),
                new_tiles_embeddings, simplify = FALSE
            ))
            colnames(new_tiles_embeddings) = colnames(tile.embeddings)
            rownames(new_tiles_embeddings) = new_tiles_metadata$id

            meta.data = rbind(meta.data, new_tiles_metadata)
            tile.shapes = c(tile.shapes, sf::st_sfc(rep(list(sf::st_multipolygon()), length(isolated_cells))))
            tile.counts = cbind(tile.counts, new_counts)
            tile.embeddings = rbind(tile.embeddings, new_tiles_embeddings)
        }
    }

    # Make Seurat object of tiles
    tile_obj = Seurat::CreateSeuratObject(
        counts = tile.counts,
        meta.data = meta.data
    )
    ## Seurat doesn't do sf shapes well
    tile_obj@meta.data$shape = tile.shapes
    row.names(tile.embeddings) = colnames(tile_obj)
    tile_obj[[reduction.name]] <- Seurat::CreateDimReducObject(embeddings = tile.embeddings, key = reduction.name)
    tile_obj[[graph.name]] = Seurat::as.Graph(res$aggs$adj)

    # Match tile_id factor levels to order in tile_obj
    obj@meta.data[[tile.id.name]] = factor(
        obj@meta.data[[tile.id.name]],
        levels=tile_obj@meta.data$id
    )

    return(list(obj=obj, tile_obj=tile_obj))
}

#' Run full DMT segmentation pipeline to make aggregated tiles from cells
#'
#' Segmentation has four main steps:
#'   1. **Preparing data structures:** A triangle mesh is constructed using Delauney
#'      triangulation and pruned to eliminate long edges. PC embeddings for each
#'      each are computed.
#'   2. **Computing gradients:** Gradients are calculated at each point by considering
#'      the difference in expression between each cell in its neighbors in the mesh.
#'      These gradients are smoothed using (anisotropic) bilateral filtering, and then
#'      gradients are defined for edges and triangles in the mesh by averaging the
#'      points that each edge or triangle contains.
#'   3. **DMT:** A scalar field is defined by taking the magnitude of the total gradient
#'      at each point/edge/triangle. Then DMT-based segmentation is performed by constructing
#'      a maximum spanning forest on the triangles and a minimum spanning forest on the points.
#'      Separatrices that separate cells into tiles of homogeneous composition are defined
#'      by tracing paths between critical points, particularly between saddle edges and maximum
#'      triangles.
#'   4. **Aggregation:** Tiles from DMT-based segmentation are merged using single-linkage
#'      agglomerative clustering to obtain tiles containing a number of cells between a
#'      user-provided minimum and maximum value. Pairs of adjacent tiles are scored according
#'      to their transcriptional similarity, compactness of shape after merging, and number
#'      of cells in order to prioritize favorable merges in each agglomerative clustering step.
#'
#' ## Computing gradients (smoothing)
#'
#' Gradient fields are smoothed using bilateral filtering,
#' in which the smoothed gradient of each point is computed as
#' the weighted average of the neighbors' gradients, considering
#' both distance in space and also similarity in gradients.
#' The weight of each neighbor is computed from the product of two scores:
#' * `distance` score: Generally, closer neighbors have greater weight.
#'   * if `'euclidean'`: Gaussian transformation of the Euclidean distance
#'     of a cell from its neighbor, so that more distant neighbors have less weight.
#'   * if `'projected'`: An anisotropic filter that accounts for expected
#'     change in expression along the direction of the neighbor. The expected
#'     change in expression is calculated from the gradient field as the total
#'     derivative in the direction of the neighbor. This change in expression is
#'     then Gaussian transformed so that neighbors that are more distant along the
#'     direction of greatest change have less weight.
#'   * if `'constant'`: All neighbors have equal `distance` weights
#' * `similarity` score: Generally, neighbors with more similar gradients have
#'   greater weight
#'   * if `'euclidean'`: Gaussian transformation of the Euclidean distance
#'     between a cell's gradient field and its neighbor's gradient field.
#'   * if `'projected'`: Gaussian transformation of the cosine distance
#'     between a cell's gradient field and its neighbor's gradient field.
#'   * if `'constant'`: All neighbors have equal `similarity` weights
#'
#' ## Aggregation (scores)
#'
#' Pairs of adjacent tiles are scored according to their transcriptional similarity,
#' compactness of shape after merging, and number of cells in order to prioritize
#' favorable merges in each agglomerative clustering step. The `dscore` for each edge
#' is computed as the product of the following three factors, where higher scores
#' favor merging of adjacent tiles:
#'   * `w`: A 2-cluster GMM is used to determine the mean `mu` and standard deviation `sig` of the distance
#'     between tiles that have similar gene expression (Euclidean distance `d` in PC space).
#'     Then we define `d_mu = mu + sig` and `d_sig = alpha * sig`, and calculate `w` as
#'     `w = 0.5 - 1 / (1 + exp(-(d - d_mu) / d_sig))`. Ranges from -0.5 to 0.5. If adjacent tiles
#'     are very dissimilar (`d >> d_mu`), then `d` is large, and `w` is close to `-0.5`. If adjacent
#'     tiles are very similar (`d < d_mu`), then `d` is small, and `w` is positive.
#'   * `score_size`: `(1 - npts_from/max_npts) * (1 - npts_to/max_npts)`. Ranges from 0 to 1.
#'     In the first round of aggregation, if merging the two tiles would have a total number
#'     of points ≥`max_npts`, then `score_size` is set to `-Inf`, which prevents merging.
#'     In the second round of aggregation, `dscore` is set to `-Inf` if both adjacent tiles
#'     have at least `min_npts` cells, to prioritize merging of small tiles.
#'   * `dC`: `.5 * (C_merge - C_from - C_to + 1)`. Ranges from 0 to 1.
#'
#' @rdname GetTiles.default
#' @family GetTiles
#'
#' @param X,Y A pair of numeric vectors with the coordinates for each of `num_cells` points.
#' @param counts A `num_genes` x `num_cells` gene-by-cell matrix of transcript counts.
#'   Optional if `embeddings` are provided directly.
#' @param embeddings A `num_cells` x `num_dim` matrix of cell embeddings across all latent dimensions.
#'   If missing, cell embeddings are calculated using PCA. If provided, the `npcs` parameter is ignored.
#' @param loadings (Optional) A `num_genes` x `num_dim` matrix of gene loadings.
#' @param meta_data A data frame with additional cell metadata to include in `dmt$pts`.
#' @param meta_vars_include Names of columns in meta_data to include in `dmt$pts`.
#' @param group.by Name of column in `meta_data` that provides the group IDs. Tessera tiles are
#'   constructed separately for each group (which could be separate experimental samples or FOVs).
#' @param npcs Number of PCs to compute for input to segmentation.
#'   Ignored if `embeddings` are provided directly.
#' @param smooth_emb Number of smoothing iterations to perform on the cell embeddings
#'   prior to gradient computation. If a vector, then embeddings after each specified
#'   iteration are concatenated. If `0` is included, then the original embeddings are also included.
#' @param prune_thresh_quantile Floating point value between 0 and 1, inclusive.
#'   Quantile of edge length above which edges are pruned. Defaults to 0.95.
#' @param prune_min_cells Minimum number of cells required for a connected
#'   component of triangles to be kept. Defaults to 10.
#' @param prune_thresh Edge length above which edges are pruned. If equal to NA,
#'   then this value is ignored and `thresh_quantile` is used to compute
#'   the threshold. Otherwise, if `thresh` is set, then `thresh_quantile`
#'   is ignored. Defaults to NA.
#' @param smooth_distance One of `c('none', 'euclidean', 'projected', 'constant')`.
#'   If either `smooth_distance` or `smooth_similarity` is `'none'`,
#'   then no smoothing of the gradient field is conducted. Defaults to `'projected'`.
#' @param smooth_similarity One of `c('none', 'euclidean', 'projected', 'constant')`.
#'   If either `smooth_distance` or `smooth_similarity` is `'none'`,
#'   then no smoothing of the gradient field is conducted. Defaults to `'projected'`.
#' @param smooth_iter Number of rounds of gradient smoothing.
#' @param max_npts Maximum number of cells allowed in each tile during the
#'   agglomerative clustering phase.
#' @param min_npts Minimum number of cells allowed in each tile during the
#'   agglomerative clustering phase.
#' @param alpha Parameter for scoring transcriptional similarity between adjacent tiles during
#'   the agglomerative clustering phase. For `alpha`, 0.2 = conservative merging, 2 = liberal merging.
#' @param future.globals.maxSize Maximum allowed size (in bytes) of global variables that are exported to each parallel worker.
#'   Increase this value if you get an error about global object size. Default is 8*1024^3 (8 GB).
#' @param consolidate Whether to consolidate results from multiple groups into a single collection of
#'   points and tiles (TRUE) or to return a list of separate results for each group (FALSE).
#' @param verbose Whether to print progress messages for each stage of the segmentation pipeline.
#'
#'
#' @returns If `consolidate==TRUE`, a List with the results of segmentation, combined across groups
#' (otherwise, if `consolidate==FALSE`, a named List, which contains separate results for each group):
#' \item{dmt}{Mesh data structures with input points/edges/triangles and the results from segmentation:
#'   * `pts`: A data table with `num_cells_pruned` rows containing cells in the mesh that
#'     remain after Delauney triangulation and pruning, with the following columns:
#'     * `X`,`Y`: Coordinates of each cell.
#'     * `ORIG_ID`: Index of each cell in the original inputs to `GetTiles()` (`X`, `Y`, `counts`, `meta_data`)
#'     * Columns from `meta_vars_include`.
#'     * `f`: Scalar value used for initial DMT-based segmentation, computed from the spatial gradient
#'       in expression at each point.
#'     * `agg_id`: Unique ID for the tile that each point belongs to in the final segmentation.
#'   * `tris`: A data table with `num_triangles` rows containing triangles in the mesh after
#'     Delauney triangulation and pruning, with the following columns:
#'     * `X`,`Y`: Coordinates of each triangle's centroid.
#'     * `area`: Area of triangle.
#'     * `height`: Largest height of each triangle.
#'     * `external`: Logical value that is `TRUE` if the triangle is a degenerate triangle
#'       that was added along a boundary edge to ensure that every edge is adjacent to two
#'       triangles. Degenerate triangles only have two vertices (the endpoints of the boundary edge).
#'     * `f`: Scalar value used for initial DMT-based segmentation, computed from the spatial gradient
#'       in expression for triangle.
#'   * `edges`: A data table with `num_edges` rows containing edges between adjacent points and triangles
#'     in the mesh after Delauney triangulation and pruning, with the following columns:
#'     * `from_pt`,`to_pt`: Indices of the two adjacent cells that are connected by an edge.
#'     * `from_tri`,`to_tri`: Indices of the two adjacent triangles that are connected by an edge.
#'     * `x0_pt`,`x1_pt`,`y0_pt`,`y1_pt`: Coordinates of the two adjacent cells that are connected by an edge.
#'     * `x0_tri`, `x1_tri`, `y0_tri`, `y1_tri`: Coordinates of the two adjacent triangles (centroids) that are connected by an edge.
#'     * `length_pt`, `length_tri`: Distance between adjacent cells or triangle centroids. Used for pruning step.
#'       (Warning - these are computed prior to pruning and are not updated after pruning and adding exterior triangles.)
#'     * `boundary`: Logical value that is `TRUE` if the edge is at the boundary and is adjacent to
#'       only a single internal triangle. For boundary edges, a degenerate external triangle is added along
#'       the boundary edge to ensure that every edge is adjacent to two triangles.
#'     * `f_prim`,`f_dual`: Scalar values used for initial DMT-based segmentation, computed from the spatial gradient
#'       in expression at each edge. Primal edges connect points and average the `f` values at the two adjacent points.
#'       Dual edges connect triangles and average the `f` values at the two adjacent triangles.
#'     * `agg_from`,`agg_to`: Unique ID for the tile that each adjacent point belongs to in the final segmentation.
#'   * `tri_to_pt`: A `num_triangles` x `num_cells_pruned` sparse matrix with value 1 at `(i,j)` if
#'     triangle `i` has point `j` as a vertex. Each internal triangle has 3 vertices, and each degenerate
#'     external triangle has 2 vertices.
#'   * `counts`: A `num_genes` x `num_cells_pruned` gene-by-cell matrix of transcript counts.
#'   * `udv_cells`: A List with the PC embeddings for each cell, which are used for segmentation.
#'     * `loadings`: If not provided as input, a `num_genes` x `npcs` matrix of gene loadings
#'       for each PC. Each column is a unit vector.
#'     * `embeddings`: If not provided as input, a `num_cells` x `npcs` matrix of cell
#'       embeddings across all PCs. Each column `j` has magnitude
#'       equal to the `j`th singular value. That is, PCs with
#'       larger contribution to the total variance will have
#'       embeddings of proportionally larger magnitude.
#'   * `prim`: The primal minimum spanning forest on points. A List with the following attributes:
#'     * `edges`: A data table with `forest_size` rows, where each row is a directed edge
#'       in the minimum spanning forest. There are six columns:
#'         * `from,to`: Index of source and target points for each edge.
#'         * `x0,y0`: Coordinates of source point for each edge.
#'         * `x1,y1`: Coordinates of target point for each edge.
#'     * `saddles`: A length `num_saddles` vector with edge indices for possible saddle edges.
#'     * `labels`: A length `num_points` vector of labels for the connected components in the
#'       minimum spanning tree. Each connected component is labeled by the index of its critical point.
#'     * `minima`: A length `num_critpts` vector of critical points (minima).
#'     * `parent`: A length `num_points` vector containing the parent (source) point for each
#'       point in the directed spanning forest. Critical points have no parent, so the value should be ignored.
#'     * `parent_edge`: A length `num_points` vector containing the directed edge that has
#'       each point as a target node. Critical points have no parent edge, so the value should be ignored.
#'   * `dual`: The dual maximum spanning forest on triangles. A List with the following attributes:
#'     * `edges`: A data.table with `forest_size` rows, where each row is a directed edge
#'       in the maximum spanning forest. There are six columns:
#'       * `from,to`: Index of source and target triangles for each edge.
#'       * `x0,y0`: Coordinates of source triangle for each edge.
#'       * `x1,y1`: Coordinates of target triangle for each edge.
#'     * `saddles`: A length `num_saddles` vector with edge indices for possible saddle edges.
#'     * `labels`: A length `num_triangles` vector of labels for the connected components in the
#'       maximum spanning tree. Each connected component is labeled by the index of its critical triangle.
#'     * `maxima`: A length `num_critpts` vector of critical triangles (maxima).
#'     * `parent`: A length `num_triangles` vector containing the parent (source) triangle for each
#'       triangle in the directed spanning forest. Critical triangles have no parent, so the value should be ignored.
#'     * `parent_edge`: A length `num_triangles` vector containing the directed edge that has
#'       each triangle as a target node. Critical triangles have no parent edge, so the value should be ignored.
#'   * `e_sep`: A length `num_sep_edges` vector of edge indices that make up the separatrices,
#'     which separate points into different components.}}
#' \item{aggs}{The tiles that result from DMT-based segmentation and agglomeration.
#'   A List data structure stores the tiles and their adjacencies using the following attributes:
#'   \itemize{
#'   * `meta_data`: A data table with `num_tiles` rows with metadata for each tile:
#'     * `ID`: Unique ID for each tile.
#'     * `X`,`Y`: Centroid of each tile.
#'     * `npts`: Number of points in each tile.
#'     * `shape`: A `sfc` list-column with the geometries for each tile.
#'     * `area,perimeter`: Area and perimeter of each tile.
#'   * `edges`: Additional attributes are calculated:
#'     * `from,to`: Tile IDs for the two tiles bordering this edge.
#'     * `x0,y0,x1,y1`: Centroid coordinates for the two tiles bordering this edge.
#'     * `area,npts`: Sum of areas and numbers of points in the two tiles bordering this edge.
#'     * `edge_length`: Total length of the border between the `from` and `to` tiles.
#'     * `dscore`: Overall score for merging two tiles. Product of `w`, `score_size`, and `dC`.
#'     * `w`: Gene expression similarity score.
#'     * `score_size`: Penalizes tiles with many points.
#'     * `perimeter_merge`: Perimeter of merged tile.
#'   * `pcs`: A `num_tiles` x `npcs` matrix with the average embedding value over all cells
#'     in each tile.
#'   * `pcs_merged`: A `num_edges` x `npcs` matrix with average PCs for the new tile if
#'     the two adjacent tiles connected by the edge were merged.
#'   * `d_mu`,`d_sig`: Parameters used to calculate `w` in the edge score `dscore`.
#'   * `aggmap`: A length `orig_num_tiles` vector mapping each original tile ID to the new
#'     tile IDs after merging.
#'   * `adj`: Sparse adjacency matrix between all tiles (if `consolidate==TRUE`).
#'   * `counts`: A `num_genes` x `num_tiles` gene-by-tile matrix of aggregated transcript counts.}
#'
#' @export
GetTiles.default = function(
    X, Y,
    counts = NULL,
    embeddings = NULL,
    loadings = NULL,
    meta_data = NULL,
    meta_vars_include = NULL,
    group.by = NULL,

    ###### STEP 0 ######
    npcs = 25,
    smooth_emb = c(0, 1),

    ## Graph pruning
    prune_thresh_quantile = 0.99,
    prune_min_cells = 1,
    prune_thresh = NA,

    ###### STEP 1: GRADIENTS ######
    smooth_distance = c('none', 'euclidean', 'projected', 'constant')[3],
    smooth_similarity = c('none', 'euclidean', 'projected', 'constant')[3],
    smooth_iter = 1,
    on_edges = FALSE,

    ###### STEP 2: DMT ######

    ###### STEP 3: AGGREGATION ######
    max_npts = 50,
    min_npts = 5,
    alpha = 1, ## 0.2 = conservative merging, 2 = liberal merging

    ## future_map parameters
    .progress = TRUE,
    .options = NULL,
    future.globals.maxSize = 8 * 1024^3,  # 8 GB

    consolidate = TRUE,
    verbose = FALSE
) {
    if (length(X) != length(Y)) {
        stop('X and Y have different lengths')
    }

    if (is.null(counts)) {
        if (is.null(embeddings)) {
            stop('Both counts and embeddings are missing. Must supply at least one.')
        }
        warning('No counts provided. Using placeholder values.')
        counts = sparseMatrix(c(), c(), dims = c(0, length(X)))
    }
    if (is.null(embeddings)) {
        warning('No embeddings provided. Calculating embeddings using PCA.')
    }

    if (is.null(group.by)) {
        warning('No value for group.by provided. Analyzing as a single sample.')
        if (is.null(meta_data)) {
            group.by = 'group'
            meta_data = data.frame(group = factor(rep(1, length(X))))
        } else {
            group.by = tail(make.unique(c(colnames(meta_data), 'group')), n = 1) # avoid overwriting existing columns
            meta_data[[group.by]] = factor(rep(1, length(X)))
        }
    } else {
        if (!(group.by %in% colnames(meta_data))) {
            stop('group.by must be a column of meta_data')
        } else if (!(group.by %in% colnames(meta_vars_include))) {
            meta_vars_include = c(meta_vars_include, group.by)
        }
    }

    # Temporarily increase global size limit for future
    old_opts <- options(future.globals.maxSize = future.globals.maxSize)
    on.exit(options(old_opts), add = TRUE)

    if (is.null(.options)) {
        .options=furrr::furrr_options(stdout = TRUE, seed = TRUE)
    }

    groups = unique(as.character(unique(meta_data[[group.by]])))
    names(groups) = groups

    res = future_map(groups, function(group) {

        idx = which(meta_data[[group.by]] == group)

        ## STEP 0: PREPARE DATA STRUCTURES
        if (verbose) message('STEP 0: PREPARE DATA STRUCTURES')
        dmt = init_data(X[idx], Y[idx], counts[,idx], meta_data[idx,], meta_vars_include)
        dmt = prune_graph(dmt, thresh_quantile = prune_thresh_quantile,
                        mincells = prune_min_cells, thresh = prune_thresh)
        dmt = add_exterior_triangles(dmt)

        if (is.null(embeddings)) {
            dmt$udv_cells = do_pca(dmt$counts, npcs)
        } else {
            dmt$udv_cells = list(
                loadings = loadings,
                embeddings = embeddings[idx,][as.integer(dmt$pts$ORIG_ID),]
            )
        }
        if (any(smooth_emb > 0)) {
            dmt = smooth_embedding(dmt, smooth_emb=smooth_emb)
        }

        ## STEP 1: GRADIENTS
        if (verbose) message('STEP 1: GRADIENTS ')
        field = compute_gradients(
            dmt, smooth_distance, smooth_similarity,
            smooth_iter = smooth_iter, on_edges = on_edges)
        field = compress_gradients_svd(field)

        ## STEP 2: DMT
        if (verbose) message('STEP 2: DMT')
        dmt = dmt_set_f(dmt, field)
        dmt$prim = do_primary_forest(dmt)
        dmt$dual = do_dual_forest(dmt)
        dmt$e_sep = dmt_get_separatrices(dmt)
        dmt = dmt_assign_tiles(dmt)
        aggs = dmt_init_tiles(dmt)

        ## STEP 3: AGGREGATION
        if (verbose) message('STEP 3: AGGREGATION')
        ## First, main aggregation
        aggs = init_scores(aggs, agg_mode=2, alpha=alpha, max_npts=max_npts)
        aggs = merge_aggs(aggs, agg_mode=2, max_npts=max_npts)
        dmt = update_dmt_aggid(dmt, aggs)
        aggs = update_agg_shapes(dmt, aggs)

        ## Then, clean up stray small aggs
        aggs = init_scores(aggs, agg_mode=3, alpha=alpha, min_npts=min_npts)
        aggs = merge_aggs(aggs, agg_mode=3, min_npts=min_npts)
        dmt = update_dmt_aggid(dmt, aggs)
        aggs = update_agg_shapes(dmt, aggs)
        aggs$meta_data$shape = st_cast(aggs$meta_data$shape, "MULTIPOLYGON")

        stopifnot(all(dmt$pts[[group.by]] == group))
        aggs$meta_data[[group.by]] = group
        dmt$pts[[group.by]] = group
        aggs$edges[[group.by]] = group

        dmt$pts$ORIG_ID = idx[dmt$pts$ORIG_ID]

        return(list(dmt=dmt, aggs=aggs))
    }, .progress=.progress, .options=.options)

    if (consolidate) {
        if (length(res) > 1) {
            res = ConsolidateResults(res, group.by)
        } else {
            res = res[[1]]
        }
        res$aggs = AddAggsAdjacencyMatrix(res$aggs)
        return(res)
    } else {
        return(res)
    }
}

#' Consolidate Tessera results from multiple samples (groups) after constructing
#' Tessera tiles separately on cells from each group.
#'
#' @param res Output of running GetTiles (when consolidate == FALSE).
#' @param group.by Name of metadata variable that identifies distinct groups.
#'
#' @export
ConsolidateResults = function(res, group.by) {
    all_aggs = list()
    all_aggs$meta_data = rbindlist(lapply(names(res), function(group) {res[[group]]$aggs$meta_data}))
    all_aggs$meta_data[[group.by]] = factor(all_aggs$meta_data[[group.by]])
    all_aggs$meta_data$id = paste0(all_aggs$meta_data[[group.by]], '_', all_aggs$meta_data$id)
    if (length(unique(all_aggs$meta_data$id)) != length(all_aggs$meta_data$id)) {
        stop('not all tile IDs are unique')
    }

    all_aggs$counts = do.call(cbind, lapply(names(res), function(group) {res[[group]]$aggs$counts}))
    colnames(all_aggs$counts) = all_aggs$meta_data$id

    all_aggs$pcs = do.call(rbind, lapply(names(res), function(group) {res[[group]]$aggs$pcs}))

    all_aggs$edges = rbindlist(lapply(names(res), function(group) {res[[group]]$aggs$edges}))
    all_aggs$edges[[group.by]] = factor(all_aggs$edges[[group.by]])
    all_aggs$edges$from = paste0(all_aggs$edges[[group.by]], '_', all_aggs$edges$from)
    all_aggs$edges$to = paste0(all_aggs$edges[[group.by]], '_', all_aggs$edges$to)

    all_dmt = list()
    all_dmt$pts = rbindlist(lapply(names(res), function(group) {res[[group]]$dmt$pts}))
    all_dmt$pts[[group.by]] = factor(all_dmt$pts[[group.by]])
    all_dmt$pts$agg_id = paste0(all_dmt$pts[[group.by]], '_', all_dmt$pts$agg_id)

    return(list(dmt=all_dmt, aggs=all_aggs))
}

#' Construct tile adjacency matrix from consolidated GetTiles output.
#'
#' @param aggs Aggregated tile information after consolidation.
#'
#' @export
AddAggsAdjacencyMatrix = function(aggs) {
    aggs$adj <- as.matrix(igraph::graph_from_data_frame(
        d = data.frame(from = aggs$edges$from, to = aggs$edges$to),
        vertices = data.frame(name = aggs$meta_data$id),
        directed = FALSE))
    stopifnot(all(colnames(aggs$adj) == aggs$meta_data$id))
    stopifnot(all(rownames(aggs$adj) == aggs$meta_data$id))
    return(aggs)
}

