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
