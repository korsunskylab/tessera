#' Log-normalization for counts data
#' 
#' Normalizes all cells to have the same total counts
#' and computes log(1+x) transform.
#' 
#' @param A A genes x cells counts matrix.
#'   Must be convertible to dgCMatrix.
#' @param scaling_factor Number of total counts to normalize
#'   all cells to. If NULL, then the median counts across
#'   all cells is used. Defaults to NULL.
#' 
#' @returns Normalized genes x cells matrix in
#'   dgCMatrix sparse matrix format.
#' 
#' @export
normalize_data = function(A, scaling_factor = NULL) {
    if (is.null(scaling_factor)) {
        scaling_factor = median(Matrix::colSums(A))
    }
    if (!'dgCMatrix' %in% class(A)) A <- as(A, "dgCMatrix")
    A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
    A@x <- scaling_factor * A@x
    A@x <- log(1 + A@x)    
	return(A)
}

#' Z-score a sparse matrix across each row or column
#' 
#' @param A A log-transformed genes x cells counts matrix.
#'   Must be convertible to dgCMatrix.
#' @param margin Subscript over which to compute Z-score.
#'   If `1`, then each row is Z-scored; otherwise, each
#'   column is Z-scored.
#' 
#' @returns Z-scored genes x cells matrix in
#'   dgCMatrix sparse matrix format.
#' 
#' @export
scale_data = function (A, margin = 1, thresh = 10) {
    A <- as(A, "dgCMatrix")
    if (margin != 1) 
        A <- t(A)
    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) 
        res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

#' Compute PCA embeddings from a raw counts matrix
#' 
#' First, log-normalizes and Z-scores the counts matrix
#' and then performs PCA using SVD.
#' 
#' @param counts A `n_genes` x `n_cells` counts matrix.
#'   Must be convertible to dgCMatrix.
#' @param npcs Number of PCs to compute.
#' 
#' @returns A list with two features:
#' * `loadings`: A `n_genes` x `npcs` matrix of gene loadings
#'   for each PC. Each column is a unit vector.
#' * `embeddings`: A `n_cells` x `npcs` matrix of cell
#'   embeddings across all PCs. Each column *j* has magnitude
#'   equal to the *j*th singular value. That is, PCs with
#'   larger contribution to the total variance will have
#'   embeddings of proportionally larger magnitude.
#' 
#' @export
do_pca = function(counts, npcs) {
    logcpx = normalize_data(counts)
    Z = scale_data(logcpx)
    Z = Z[which(is.na(Matrix::rowSums(Z)) == 0), ]
    pres = RSpectra::svds(Z, npcs)
    V = sweep(pres$v, 2, pres$d, '*')
    colnames(V) = paste0('PC', 1:npcs)
    row.names(V) = colnames(Z)
    colnames(pres$u) = paste0('PC', 1:npcs)
    row.names(pres$u) = row.names(Z)
    return(list(loadings = pres$u, embeddings = V))
}

#' Smooth embeddings along edges
#' 
#' @param dmt A DMT object.
#' @param smooth_emb Number of smoothing iterations to perform.
#'   If a vector, then embeddings after each specified
#'   iteration are concatenated. If `0` is included, then
#'   the original embeddings are also included.
#'
#' @returns The input DMT object with smoothed embeddings stored in
#'   `dmt$udv_cells$embeddings`.
#' 
#' @export
smooth_embedding = function(
    dmt,
    smooth_emb = 0
) {

    adj = as.matrix(igraph::from_edgelist()$fun(as.matrix(dmt$edges[, .(from_pt, to_pt)]), directed=FALSE))
    diag(adj) = 1
    adj = adj / Matrix::colSums(adj)  # normalize

    embeddings = dmt$udv_cells$embeddings

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
    rownames(embeddings) = rownames(dmt$udv_cells$embeddings)
    colnames(embeddings) = paste0('PC_', 1:ncol(embeddings))
    dmt$udv_cells$embeddings = embeddings

    if (!is.null(dmt$udv_cells$loadings)) {
        loadings = dmt$udv_cells$loadings

        loadings = do.call(cbind, replicate(length(smooth_emb), loadings, simplify=FALSE))
        rownames(loadings) = rownames(dmt$udv_cells$loadings)
        colnames(loadings) = paste0('PC_', 1:ncol(embeddings))

        dmt$udv_cells$loadings = loadings
    }

    return(dmt)
}
