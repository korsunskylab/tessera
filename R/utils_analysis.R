#' Run UMAP and save fgraph and embeddings in Seurat object
#' 
#' @param obj A Seurat object.
#' @param reduction Name of dimensional reduction in `obj` to use as input to UMAP.
#' @param dims Dimensions of `reduction` to use as input to UMAP. If NULL, use all dimensions.
#' @param fgraph_only If TRUE, only compute and store the fuzzy simplicial graph (fgraph) and skip UMAP embedding.
#' @param graph.name Name of graph to store the fgraph in `obj`. Defaults to `<assay>_fgraph`.
#' @param reduction.name Name of dimensional reduction to store the UMAP embeddings in `obj`. Defaults to "umap".
#' @param assay Assay to set as default assay for the new dimensional reduction. If NULL, use the default assay of `obj`.
#' @param key Key prefix to use for the new dimensional reduction. Defaults to "UMAP_".
#' @param n_neighbors Number of nearest neighbors to use in UMAP. See `?uwot::umap` for details.
#' @param n_components Number of UMAP dimensions to compute. Ignored if `fgraph_only` is TRUE.
#' @param metric Distance metric to use in UMAP. See `?uwot::umap` for details.
#' @param spread UMAP spread parameter. See `?uwot::umap` for details.
#' @param min_dist UMAP minimum distance parameter. See `?uwot::umap` for details.
#' @param n_threads Number of threads to use in UMAP. See `?uwot::umap` for details.
#' @param fast_sgd Whether to use the fast stochastic gradient descent optimization in UMAP. See `?uwot::umap` for details.
#' @param verbose Whether to print progress messages.
#' @param ... Additional parameters to pass to `uwot::umap`.
#' 
#' @returns The input Seurat object with the fgraph and (optionally) UMAP embeddings added.
#'
#' @export
RunUMAPCustom = function(
    obj, reduction = "pca", dims = NULL, fgraph_only = FALSE,
    graph.name = NULL, reduction.name = "umap", assay = NULL, key = "UMAP_",
    n_neighbors = 30, n_components = 2, metric = "cosine", spread = 1, min_dist = 0.3,
    n_threads = NULL, fast_sgd = TRUE, verbose = TRUE, ...
) {
    if (is.null(assay)) {
        assay = Seurat::DefaultAssay(obj)
    }
    if (is.null(graph.name)) {
        graph.name = paste0(assay, '_fgraph')
    }

    embeddings = Seurat::Embeddings(obj, reduction)
    if (!is.null(dims)) {
        embeddings = embeddings[,dims]
    }

    if (fgraph_only) {
        out = list()
        out$fgraph = uwot::similarity_graph(
            X = embeddings,
            n_neighbors = n_neighbors,
            metric = metric,
            method = "umap",
            verbose = verbose
        )
    } else {
        out = uwot::umap(
            X=embeddings, ret_extra = "fgraph",
            n_threads = n_threads, fast_sgd = fast_sgd, verbose = verbose,
            n_neighbors = n_neighbors, n_components = n_components,
            metric = metric, spread = spread, min_dist = min_dist, ...
        )
        colnames(out$embedding) = paste0("UMAP_", seq_len(n_components))
        obj[[reduction.name]] = Seurat::CreateDimReducObject(embeddings = out$embedding,
                                                     key = key, assay = assay)
    }
    rownames(out$fgraph) = colnames(obj)
    colnames(out$fgraph) = colnames(obj)
    obj[[graph.name]] = Seurat::as.Graph(out$fgraph)

    return(obj)
}

#' Merge clusters in a Seurat object
#' 
#' @param obj A Seurat object.
#' @param to_merge A vector of cluster labels to merge.
#' @param new_label The new label for the merged cluster.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels
#'  to be merged. Defaults to 'seurat_clusters'.
#' @param new.clusters.name Name of the metadata column to store the new cluster labels.
#'  If NULL, defaults to `<clusters.name>_m<new_label>`.
#' 
#' @returns The input Seurat object with the merged cluster labels added to metadata.
#' 
#' @export
MergeClusters = function(obj, to_merge, new_label, clusters.name = 'seurat_clusters', new.clusters.name = NULL) {
    stopifnot(is.factor(obj@meta.data[[clusters.name]]))
    
    if (is.null(new.clusters.name)) {
        new.clusters.name = paste0(clusters.name, '_m', new_label)
    }

    old_levels = levels(obj@meta.data[[clusters.name]])
    new_levels = c(old_levels[!(old_levels %in% to_merge)], new_label)

    obj@meta.data[[new.clusters.name]] = as.character(obj@meta.data[[clusters.name]])
    obj@meta.data[[new.clusters.name]][obj@meta.data[[new.clusters.name]] %in% to_merge] = new_label

    obj@meta.data[[new.clusters.name]] = factor(
        obj@meta.data[[new.clusters.name]],
        levels = new_levels
    )

    return(obj)
}

#' Add sub-cluster labels to a Seurat object
#' 
#' @param obj A Seurat object.
#' @param obj_sub A Seurat object containing the sub-clustered cells.
#' @param cluster The parent cluster label in `obj` that was sub-clustered.
#' @param clusters.name Name of the metadata column in `obj` containing the parent cluster labels.
#' @param sub.clusters.name Name of the metadata column to store the sub-cluster labels.
#'  If NULL, defaults to `<clusters.name>_s<cluster>`.
#'
#' @returns The input Seurat object with sub-cluster labels added to metadata.
#' 
#' @export
AddSubClusterLabels = function(
    obj, obj_sub, cluster, clusters.name = "seurat_clusters", sub.clusters.name = NULL
) {
    stopifnot(is.factor(obj@meta.data[[clusters.name]]))
    stopifnot(sum(obj@meta.data[[clusters.name]] == cluster) == ncol(obj_sub))
    stopifnot(all(rownames(obj@meta.data %>%
                           select(!!sym(clusters.name)) %>%
                           filter(!!sym(clusters.name) == cluster)) ==
                  rownames(obj_sub@meta.data)))

    if (is.null(sub.clusters.name)) {
        sub.clusters.name = paste0(clusters.name, '_s', cluster)
    }
    
    clusters.old = obj@meta.data[[clusters.name]]
    obj@meta.data[[sub.clusters.name]] = as.character(clusters.old)
    obj@meta.data[rownames(obj_sub@meta.data), sub.clusters.name] = paste0(cluster, '_', obj_sub$seurat_clusters)
    idx = which(levels(clusters.old) == cluster)
    
    new_levels = c(
        levels(clusters.old)[seq_len(idx-1)],
        paste0(cluster, '_', levels(obj_sub$seurat_clusters)),
        levels(clusters.old)[idx + seq_len(nlevels(clusters.old)-idx)]
    )
    obj@meta.data[[sub.clusters.name]] = factor(obj@meta.data[[sub.clusters.name]], levels = new_levels)

    return(obj)
}

#' Create a Seurat object for sub-clustering a specific cluster
#' 
#' @param obj A Seurat object.
#' @param cluster The cluster label to sub-cluster.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels.
#' @param npcs Number of principal components to use for sub-clustering.
#' @param fast_sgd Whether to use fast SGD in UMAP.
#' @param n_neighbors Number of neighbors to use in UMAP.
#' @param scale.factor Scale factor for normalization. If NULL, use median nCount_RNA.
#' @param use.existing.embeddings Name of existing dimensional reduction in `obj`
#'   to use for sub-clustering. If NULL, compute PCA on the subsetted data.
#' @param meta.vars.include Metadata variables to include in the sub-clustered object.
#' @param harmony.group.by.vars Metadata variables to use for Harmony integration.
#' @param early_stop Whether to use early stopping in Harmony.
#' @param ... Additional parameters to pass to `harmony::RunHarmony`. Ignored if `harmony.group.by.vars` is NULL.
#' 
#' @returns A Seurat object for the sub-clustered cells.
#'
#' @export
MakeSubClusterObj = function(
    obj, cluster,
    clusters.name = "seurat_clusters", npcs=30,
    fast_sgd = TRUE,
    n_neighbors = 15,
    assay = NULL,
    scale.factor = NULL,
    use.existing.embeddings = NULL,
    meta.vars.include = NULL,
    harmony.group.by.vars = NULL,
    early_stop = TRUE,
    ...
) {
    if (is.null(assay)) {
        assay = Seurat::DefaultAssay(obj)
    }
    
    subset_idx = which(obj@meta.data[[clusters.name]] == cluster)
    meta.vars.include = unique(c(meta.vars.include, harmony.group.by.vars))
    # counts = obj[["RNA"]]$counts[,subset_idx]
    # colnames(counts) = NULL
    obj_sub = Seurat::CreateSeuratObject(
        counts = obj[[assay]]$counts[,subset_idx],
        # counts = counts,
        meta.data = obj@meta.data[subset_idx, meta.vars.include]
    )

    if (is.null(use.existing.embeddings)) {
        if (is.null(scale.factor)) {
            scale.factor = median(obj_sub$nCount_RNA)
        }
        obj_sub <- Seurat::NormalizeData(object = obj_sub, scale.factor = scale.factor)
        obj_sub <- Seurat::FindVariableFeatures(object = obj_sub)
        obj_sub <- Seurat::ScaleData(object = obj_sub)  # might get an error in this line if the cells don't have good names
        obj_sub <- Seurat::RunPCA(object = obj_sub, npcs = npcs)
        reduction = "pca"
    } else {
        obj_sub[[use.existing.embeddings]] <- Seurat::CreateDimReducObject(
            embeddings = Seurat::Embeddings(obj, use.existing.embeddings)[subset_idx,],
            loadings = Seurat::Loadings(obj, use.existing.embeddings),
            key = use.existing.embeddings
        )
        reduction = use.existing.embeddings
    }

    if (!is.null(harmony.group.by.vars)) {
        obj_sub = harmony::RunHarmony(
            obj_sub, harmony.group.by.vars,
            early_stop = early_stop,
            reduction.use = reduction,
            reduction.save = paste0(reduction, '_harmony'),
            ...
        )
        reduction = paste0(reduction, '_harmony')
    }

    obj_sub = RunUMAPCustom(obj_sub, reduction = reduction,
                            n_neighbors = n_neighbors, fast_sgd = fast_sgd)
    return(obj_sub)
}

#' Create a Seurat object for sub-clustering a specific cluster from matching tile and cell data
#' 
#' @param tile_obj A Seurat object containing tile-level data.
#' @param cluster The cluster label to sub-cluster.
#' @param cell_obj A Seurat object containing cell-level data with a `tile_id` metadata column.
#' @param clusters.name Name of the metadata column in `tile_obj` containing the cluster labels.
#' @param npcs Number of principal components to use for sub-clustering.
#' @param fast_sgd Whether to use fast SGD in UMAP.
#' @param n_neighbors Number of neighbors to use in UMAP.
#' @param tile_assay Assay in `tile_obj` to use for sub-clustering. If NULL, use DefaultAssay.
#' @param cell_assay Assay in `cell_obj` to use for computing cell embeddings. If NULL, use DefaultAssay.
#' @param scale.factor Scale factor for normalization. If NULL, use median nCount_RNA.
#' @param use.existing.embeddings Name of existing dimensional reduction in `tile_obj`
#'   to use for sub-clustering. If NULL, compute PCA on the subsetted cell-level data.
#' @param smooth_emb Number of smoothing iterations to perform on the cell embeddings,
#'   as for `GetTiles`. If a vector, then embeddings after each specified
#'   iteration are concatenated. If `0` is included, then the original embeddings are also included.
#' @param graph.name.cells Name of the graph in `cell_obj` containing the cell adjacency graph.
#' @param meta.vars.include Metadata variables to include in the sub-clustered object.
#' @param harmony.group.by.vars Metadata variables to use for Harmony integration.
#' @param early_stop Whether to use early stopping in Harmony.
#'
#' @returns A Seurat object for the tile subset that can be used for sub-clustering.
#' 
#' @export
MakeTileSubClusterObj = function(
    tile_obj, cluster, cell_obj,
    clusters.name = "seurat_clusters", npcs=30,
    tile_assay = NULL, cell_assay = NULL,
    fast_sgd = TRUE,
    n_neighbors = 15,
    scale.factor = NULL,
    use.existing.embeddings = NULL,
    smooth_emb = c(0, 1),
    graph.name.cells = 'cell_adj',
    meta.vars.include = NULL,
    harmony.group.by.vars = NULL,
    early_stop = TRUE
) {
    stopifnot(all(levels(cell_obj$tile_id) == colnames(tile_obj)))
    stopifnot('tile_id' %in% colnames(cell_obj@meta.data))
    
    if (is.null(tile_assay)) {
        tile_assay = Seurat::DefaultAssay(tile_obj)
    }
    
    subset_idx = which(tile_obj@meta.data[[clusters.name]] == cluster)
    meta.vars.include = unique(c(meta.vars.include, harmony.group.by.vars))
    tile_obj_sub = Seurat::CreateSeuratObject(
        counts = tile_obj[[tile_assay]]$counts[,subset_idx],
        meta.data = tile_obj@meta.data[subset_idx, meta.vars.include, drop=FALSE]
    )
    
    if (is.null(use.existing.embeddings)) {
        if (is.null(cell_assay)) {
            cell_assay = Seurat::DefaultAssay(cell_obj)
        }
        
        cell_subset_idx = which(cell_obj$tile_id %in% colnames(tile_obj_sub))
        meta.vars.include = unique(c(meta.vars.include, harmony.group.by.vars, 'tile_id'))
        
        cell_obj_sub = Seurat::CreateSeuratObject(
            counts = cell_obj[[cell_assay]]$counts[,cell_subset_idx],
            meta.data = cell_obj@meta.data[cell_subset_idx, meta.vars.include, drop=FALSE]
        )
        
        if (is.null(scale.factor)) {
            scale.factor = median(cell_obj_sub$nCount_RNA)
        }
        cell_obj_sub <- Seurat::NormalizeData(object = cell_obj_sub, scale.factor = scale.factor)
        cell_obj_sub <- Seurat::FindVariableFeatures(object = cell_obj_sub)
        cell_obj_sub <- Seurat::ScaleData(object = cell_obj_sub)  # might get an error in this line if the cells don't have good names
        cell_obj_sub <- Seurat::RunPCA(object = cell_obj_sub, npcs = npcs)
        reduction = "pca"
    
        if (!is.null(harmony.group.by.vars)) {
            cell_obj_sub = harmony::RunHarmony(
                cell_obj_sub, harmony.group.by.vars,
                early_stop = early_stop,
                reduction.use = reduction,
                reduction.save = paste0(reduction, '_harmony'))
            reduction = paste0(reduction, '_harmony')
        }

        # smooth embeddings
        embeddings = Seurat::Embeddings(cell_obj_sub, reduction)
        loadings = Seurat::Loadings(cell_obj_sub, reduction)
        if (!all(smooth_emb == 0)) {
            adj = Seurat::as.sparse(cell_obj[[graph.name.cells]])[cell_subset_idx,cell_subset_idx]
            diag(adj) = 1
            adj = adj / Matrix::colSums(adj)  # normalize

            smoothed_embeddings = list()
            if (0 %in% smooth_emb) {
                smoothed_embeddings[['0']] = as.matrix(embeddings)
            }
            for (i in seq_len(max(smooth_emb))) {
                embeddings = adj %*% embeddings
                if (i %in% smooth_emb) {
                    smoothed_embeddings[[as.character(i)]] = as.matrix(embeddings)
                }
            }

            embeddings = do.call(cbind, smoothed_embeddings)
            rownames(embeddings) = colnames(cell_obj_sub)
            colnames(embeddings) = paste0('PC_', 1:ncol(embeddings))

            loadings = do.call(cbind, replicate(length(smooth_emb), loadings, simplify=FALSE))
            rownames(loadings) = rownames(Seurat::Loadings(cell_obj_sub, reduction))
            colnames(loadings) = paste0('PC_', 1:ncol(embeddings))
        }
        
        tile_obj_sub[[reduction]] <- Seurat::CreateDimReducObject(
            embeddings = aggregate_embeddings(
                embeddings, droplevels(cell_obj_sub$tile_id)
            )[colnames(tile_obj_sub),],
            loadings = loadings,
            key = reduction
        )
    } else {
        tile_obj_sub[[use.existing.embeddings]] <- Seurat::CreateDimReducObject(
            embeddings = Seurat::Embeddings(tile_obj, use.existing.embeddings)[subset_idx,],
            loadings = Seurat::Loadings(tile_obj, use.existing.embeddings),
            key = use.existing.embeddings
        )
        reduction = use.existing.embeddings
    }
    
    tile_obj_sub = RunUMAPCustom(
        tile_obj_sub, reduction = reduction,
        n_neighbors = n_neighbors, fast_sgd = fast_sgd
    )
    return(tile_obj_sub)
}

#' Find sub-clusters within a specific cluster of a Seurat object
#' 
#' @param obj A Seurat object.
#' @param cluster The cluster label to sub-cluster.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels.
#' @param sub.clusters.name Name of the metadata column to store the sub-cluster labels.
#'  If NULL, defaults to `<clusters.name>_s<cluster>`.
#' @param resolution Resolution parameter for clustering.
#' @param algorithm Clustering algorithm to use. See `?Seurat::FindClusters` for details.
#' @param npcs Number of principal components to use for sub-clustering.
#' @param method Method to use for clustering. See `?Seurat::FindClusters` for details.
#' @param n_neighbors Number of neighbors to use in UMAP.
#' @param fast_sgd Whether to use fast SGD in UMAP.
#' @param scale.factor Scale factor for normalization. If NULL, use median nCount_RNA.
#' @param use.existing.embeddings Name of existing dimensional reduction in `obj`
#'   to use for sub-clustering. If NULL, compute PCA on the subsetted data.
#' @param meta.vars.include Metadata variables to include in the sub-clustered object.
#' @param harmony.group.by.vars Metadata variables to use for Harmony integration.
#' @param early_stop Whether to use early stopping in Harmony.
#' @param return_obj_sub If TRUE, return a list with the updated `obj` and the sub-clustered object.
#' @param ... Additional parameters to pass to `harmony::RunHarmony`. Ignored if `harmony.group.by.vars` is NULL.
#' 
#' @returns The input Seurat object with sub-cluster labels added to metadata.
#' 
#' @export
FindSubClusterCustom = function(
    obj, cluster,
    clusters.name = "seurat_clusters",
    sub.clusters.name = NULL,
    resolution = 0.5, algorithm = 1, npcs=30,
    method = 'igraph',
    n_neighbors = 15,
    fast_sgd = TRUE,
    scale.factor = NULL,
    use.existing.embeddings = NULL,
    meta.vars.include = NULL,
    harmony.group.by.vars = NULL,
    early_stop = TRUE,
    return_obj_sub = FALSE,
    ...
) {
    obj_sub = MakeSubClusterObj(obj, cluster, clusters.name = clusters.name,
                                npcs = npcs, n_neighbors = n_neighbors,
                                scale.factor = scale.factor,
                                use.existing.embeddings = use.existing.embeddings,
                                meta.vars.include = meta.vars.include,
                                harmony.group.by.vars = harmony.group.by.vars,
                                early_stop = early_stop,
                                fast_sgd = fast_sgd,
                                ...
    )
    if (algorithm == 4) {
        suppressWarnings({    # the igraph conversion spits out a lot of warnings, which is slow if printed
            obj_sub <- Seurat::FindClusters(object = obj_sub, resolution = resolution, algorithm = algorithm,
                            method = method, graph.name = 'RNA_fgraph')
        })
    } else {
        obj_sub <- Seurat::FindClusters(object = obj_sub, resolution = resolution, algorithm = algorithm,
                            method = method, graph.name = 'RNA_fgraph')
    }

    obj = AddSubClusterLabels(obj, obj_sub, cluster, clusters.name = clusters.name, sub.clusters.name = sub.clusters.name)

    if (return_obj_sub) {
        list('obj' = obj, 'obj_sub' = obj_sub)
    } else {
        return(obj)
    }
}

#' Find sub-clusters within a specific cluster of a Seurat object using matching tile and cell data
#' 
#' @param obj A Seurat object containing tile-level data.
#' @param cluster The cluster label to sub-cluster.
#' @param cell_obj (Opional) A Seurat object containing cell-level data with a `tile_id` metadata column.
#'   If NULL, sub-clustering is performed using only `obj` without recomputing tile embeddings from cell data.
#' @param clusters.name Name of the metadata column in `obj` containing the cluster labels.
#' @param tile_assay Assay in `obj` to use for sub-clustering. If NULL, use DefaultAssay.
#' @param cell_assay Assay in `cell_obj` to use for computing cell embeddings. If NULL, use DefaultAssay.
#' @param sub.clusters.name Name of the metadata column to store the sub-cluster labels.
#'  If NULL, defaults to `<clusters.name>_s<cluster>`.
#' @param resolution Resolution parameter for clustering.
#' @param algorithm Clustering algorithm to use. See `?Seurat::FindClusters` for details.
#' @param npcs Number of principal components to use for sub-clustering.
#' @param method Method to use for clustering. See `?Seurat::FindClusters` for details.
#' @param n_neighbors Number of neighbors to use in UMAP.
#' @param fast_sgd Whether to use fast SGD in UMAP.
#' @param scale.factor Scale factor for normalization. If NULL, use median nCount_RNA.
#' @param use.existing.embeddings Name of existing dimensional reduction in `obj`
#'   to use for sub-clustering. If NULL, compute PCA on the subsetted cell-level data.
#' @param smooth_emb Number of smoothing iterations to perform on the cell embeddings,
#'   as for `GetTiles`. If a vector, then embeddings after each specified
#'   iteration are concatenated. If `0` is included, then the original embeddings are also included.
#' @param graph.name.cells Name of the graph in `cell_obj` containing the cell adjacency graph.
#' @param meta.vars.include Metadata variables to include in the sub-clustered object.
#' @param harmony.group.by.vars Metadata variables to use for Harmony integration.
#' @param early_stop Whether to use early stopping in Harmony.
#' @param return_obj_sub If TRUE, return a list with the updated `obj` and the sub-clustered object.
#'
#' @returns The input Seurat object with sub-cluster labels added to metadata.
#' 
#' @export
FindTileSubCluster = function(
    obj, cluster, cell_obj = NULL,
    clusters.name = "seurat_clusters",
    tile_assay = NULL, cell_assay = NULL,
    sub.clusters.name = NULL,
    resolution = 0.2, algorithm = 4, npcs=20,
    method = 'igraph',
    n_neighbors = 25,
    fast_sgd = TRUE,
    scale.factor = NULL,
    use.existing.embeddings = NULL,
    smooth_emb = c(0, 1),
    graph.name.cells = 'cell_adj',
    meta.vars.include = NULL,
    harmony.group.by.vars = NULL,
    early_stop = TRUE,
    return_obj_sub = FALSE
) {
    if (is.null(cell_obj)) {
        obj_sub = MakeSubClusterObj(
            obj, cluster, clusters.name = clusters.name,
            npcs = npcs, n_neighbors = n_neighbors,
            assay = tile_assay,
            scale.factor = scale.factor,
            use.existing.embeddings = use.existing.embeddings,
            meta.vars.include = meta.vars.include,
            harmony.group.by.vars = harmony.group.by.vars,
            early_stop = early_stop,
            fast_sgd = fast_sgd
        )
    } else {
        obj_sub = MakeTileSubClusterObj(
            obj, cluster, cell_obj,
            clusters.name = clusters.name,
            tile_assay = tile_assay, cell_assay = cell_assay,
            npcs = npcs, n_neighbors = n_neighbors,
            scale.factor = scale.factor,
            use.existing.embeddings = use.existing.embeddings,
            smooth_emb = smooth_emb,
            graph.name.cells = graph.name.cells,
            meta.vars.include = meta.vars.include,
            harmony.group.by.vars = harmony.group.by.vars,
            early_stop = early_stop,
            fast_sgd = fast_sgd
        )
    }
    
    if (algorithm == 4) {
        suppressWarnings({    # the igraph conversion spits out a lot of warnings, which is slow if printed
            obj_sub <- Seurat::FindClusters(object = obj_sub, resolution = resolution, algorithm = algorithm,
                            method = method, graph.name = 'RNA_fgraph')
        })
    } else {
        obj_sub <- Seurat::FindClusters(object = obj_sub, resolution = resolution, algorithm = algorithm,
                            method = method, graph.name = 'RNA_fgraph')
    }

    obj = AddSubClusterLabels(obj, obj_sub, cluster, clusters.name = clusters.name, sub.clusters.name = sub.clusters.name)

    if (return_obj_sub) {
        list('obj' = obj, 'obj_sub' = obj_sub)
    } else {
        return(obj)
    }
}

#' @param labels A factor of length num_cells (can be named)
#'   assigning each cell to some group (tile ID, cell type, etc.).
#'   Cells that have NA value for the group are not assigned to any group.
#' 
#' @returns A `num_groups` x `num_cells` sparse matrix with one-hot assignment
#'   of each cell to a group (unassigned if NA). Rownames are the factor
#'   levels (group names), and colnames are the factor names (cell names).
#'
#' @export
group_matrix = function(groups) {
    stopifnot(is.factor(groups))

    groups_as_int = as.integer(groups)
    not_na_idx = which(!is.na(groups_as_int))
    groups_as_int = groups_as_int[not_na_idx]

    out = Matrix::sparseMatrix(
        i = groups_as_int,
        j = not_na_idx,
        x = 1,
        dims = c(nlevels(groups), length(groups)),
        dimnames = list(levels(groups), names(groups))
    )
    return(out)
}

#' @param embeddings A `num_cells` x `embedding_dim` (or 
#'   `embeddings_dim` x `num_cells` if transposed is TRUE) matrix
#'   of cell embeddings (e.g. transcript counts, PCs, harmonized PCs,
#'   scVI latent dimensions, etc.)
#' @param groups A factor of length num_cells (can be named)
#'   assigning each cell to some category (tile ID, cell type, etc.)
#' @param mean If `TRUE`, then aggregated by computing the mean
#'   of each group. Otherwise, computes the sum.
#' @param as_matrix If `TRUE`, then casts output as dense matrix
#' @param transposed If TRUE, then `embeddings` input and output
#'   have dimensions `embeddings_dim` x `num_cells` and
#'   `embedding_dim` x `num_groups`
#'
#' @returns A `num_groups` x `embedding_dim` matrix (or 
#'   `embedding_dim` x `num_groups` if tranposed is TRUE)'
#' 
#' @export
aggregate_embeddings = function(embeddings, groups, mean = TRUE, as_matrix = TRUE, transposed = FALSE) {
    
    group_by_cell = group_matrix(groups)
    if (mean) {
        group_counts = Matrix::rowSums(group_by_cell)
        if (any(group_counts == 0)) {
            message('WARNING: Some groups have no members, will have NA embeddings')
        }
        # group_counts[group_counts == 0] = 1  # avoid divide by zero
        group_by_cell = group_by_cell / group_counts
    }
    if (transposed) {
        aggregated = Matrix::tcrossprod(embeddings, group_by_cell)
    } else {
        aggregated = group_by_cell %*% embeddings    # sum of embeddings
    }
    if (as_matrix) {
        aggregated = as.matrix(aggregated)
    }
    return(aggregated)
}

#' Add tile-level metadata to a matching cell-level Seurat object
#' 
#' @param obj A Seurat object containing cell-level data with a `tile_id` metadata column.
#' @param tile_obj A Seurat object containing tile-level data.
#' @param ... Tidy expressions specifying which metadata columns from `tile_obj`
#'   to add to `obj`. For example, `tile_size = npts`
#' @param by A named character vector specifying how to join `obj` and `tile_obj`.
#'   Defaults to `c('tile_id' = 'id')`, where `tile_id` is the metadata column
#'   in `obj` and `id` is the metadata column in `tile_obj`.
#'
#' @returns The input Seurat object with additional tile-level metadata added to `@meta.data`.
#'
#' @export
AddTileMetadata = function(
    obj, tile_obj, ...,
    by = c('tile_id' = 'id')
) {
    # capture expressions
    exprs <- enquos(...)

    if (length(exprs) == 0) {
        stop("Please supply one or more tidy expressions, e.g. tile_size = npts")
    }

    obj@meta.data = cbind(
        obj@meta.data,
        obj@meta.data %>% select(tile_id) %>% left_join(
            tile_obj@meta.data %>% select(id, !!!exprs),
            by = by
        ) %>% select(-tile_id)
    )

    return(obj)
}
