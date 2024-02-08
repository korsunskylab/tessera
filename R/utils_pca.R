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
