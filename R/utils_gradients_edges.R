#' Bilateral / anisotropic filtering of gradient field
#'
#' @export
smooth_field_edges = function(
    edges, field, edges_svd, coords, adj_idx,
    distance="projected", similarity="euclidean"
) {

    ## method to compute distance between points
    if (distance == "euclidean") {
        ## Euclidean distance
        dist_mode = 0
    } else if (distance == "projected") {
        ## Projected distance (anisotropic)
        dist_mode = 1
    } else if (distance == "constant") {
        dist_mode = 2
    } else {
        stop("Invalid distance method")
    }

    ## compute similarity between fields
    if (similarity == "euclidean") {
        ## Euclidean distance
        sim_mode = 0
    } else if (similarity == "projected") {
        sim_mode = 1
    } else if (similarity == "constant") {
        sim_mode = 2
    } else {
        stop("Invalid similarity method")
    }

    res = smooth_field_edges_cpp(
        edges$from_pt-1, # input is 0-indexed
        edges$to_pt-1, # input is 0-indexed
        field,   # arma::cube shape (2, D, E)
        edges_svd, # u, s, v
        coords,   # arma::mat shape (N, 2)
        adj_idx,    # dgCMatrix shape (N, N): adj_idx@x should store mapping
                    # to 1-indexed edge idx
        dist_mode,
        sim_mode
    )

    return(res)
}

#' Compute spatial gradient field for input to DMT
#'
#' @export
compute_gradients_edges = function(
    dmt,
    smooth_distance='projected',
    smooth_similarity='euclidean',
    smooth_iter=3
) {
    field = list()

    ## First on edges
    coords = as.matrix(dmt$pts[, .(X, Y)])
    adj_idx = Matrix::sparseMatrix(
        i = c(dmt$edges$from_pt, dmt$edges$to_pt),
        j = c(dmt$edges$to_pt, dmt$edges$from_pt),
        x = c(seq_len(nrow(dmt$edges)),
              seq_len(nrow(dmt$edges))),  # 1-indexed edges
        dims = c(nrow(dmt$pts), nrow(dmt$pts))
    )
    stopifnot(all(Matrix::diag(adj_idx) == 0))

    embeddings = dmt$udv_cells$embeddings
    field$edges = estimate_field_edges_cpp(
        coords,
        dmt$udv_cells$embeddings,
        dmt$edges$from_pt-1,
        dmt$edges$to_pt-1
    )
    field$edges_unsmoothed = field$edges
    field$edges_svd = svd_field_cpp(field$edges)

    ## Smooth edge gradient first
    for (i in seq_len(smooth_iter)) {
        field$edges = smooth_field_edges(
            edges = dmt$edges,
            field = field$edges,
            edges_svd = field$edges_svd,
            coords = coords,
            adj_idx = adj_idx,
            distance=smooth_distance,
            similarity=smooth_similarity
        )
        field$edges_svd = svd_field_cpp(field$edges)
    }

    ## Field on tris: average of gradients for boundary edges
    tri_to_edge = Matrix::sparseMatrix(
        i = c(dmt$edges$from_tri, dmt$edges$to_tri),
        j = c(seq_len(nrow(dmt$edges)), seq_len(nrow(dmt$edges))),
        x = 1,
        dims = c(nrow(dmt$tris), nrow(dmt$edges))
    )  # every edge has 2 triangles, every triangle has 3 edges, or 1 if degenerate
    tri_to_edge = tri_to_edge / Matrix::rowSums(tri_to_edge)

    field$tris = array(dim = c(2, dim(field$edges)[2], nrow(dmt$tris)))
    field$tris[1,,] = as.matrix(Matrix::tcrossprod(field$edges[1,,], tri_to_edge))
    field$tris[2,,] = as.matrix(Matrix::tcrossprod(field$edges[2,,], tri_to_edge))

    ## Field on pts: average of gradients for incoming edges
    pt_to_edge = Matrix::sparseMatrix(
        i = c(dmt$edges$from_pt, dmt$edges$to_pt),
        j = c(seq_len(nrow(dmt$edges)), seq_len(nrow(dmt$edges))),
        x = 1,
        dims = c(nrow(dmt$pts), nrow(dmt$edges))
    )  # every edge has 2 pts, each pt can have a range of neighbors
    pt_to_edge = pt_to_edge / Matrix::rowSums(pt_to_edge)

    field$pts = array(dim = c(2, dim(field$edges)[2], nrow(dmt$pts)))
    field$pts[1,,] = as.matrix(Matrix::tcrossprod(field$edges[1,,], pt_to_edge))
    field$pts[2,,] = as.matrix(Matrix::tcrossprod(field$edges[2,,], pt_to_edge))

    ## Backwards compatibility
    field$edges_pts = field$edges
    field$edges_tris = field$edges

    return(field)
}
