#' Bilateral / anisotropic filtering of gradient field
#'
#' Gradient fields are smoothed using bilateral filtering,
#' in which the smoothed gradient of each point is computed as
#' the weighted average of the neighbors' gradients, considering
#' both distance in space and also similarity in gradients.
#'
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
#' @param coords A `N` x `2` matrix of cell coordinates.
#' @param field A `2` x `D` x `N` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` latent variables at every point in space.
#' @param adj A `N` x `N` sparse adjacency matrix
#'   in dgCMatrix format.
#' @param include_self A boolean whether or not to include the
#'   each point's gradient in its own smoothed value. Defaults to TRUE.
#' @param distance Method for computing distance score in weighted average.
#'   See description for details. Defaults to `'euclidean'`.
#' @param similarity Method for computing similarity score in weighted average.
#'   See description for details. Defaults to `'euclidean'`.
#'
#' @returns A `2` x `D` x `N` array in column-major ordering
#'   containing the smoothed spatial gradient in expression for each of
#'   `D` latent variables at every point in space.
#'
#' @export
smooth_field = function(coords, field, adj, include_self=TRUE,
                        distance="euclidean", similarity="euclidean") {
    if (include_self) {
        diag(adj) = 1
    } else {
        diag(adj) = 0
    }

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

    res = smooth_field_cpp(
        pvec = diff(adj@p),
        adj_i = adj@i, ## keep it 0-indexed
        adj_p = adj@p,
        field = field,
        coords = coords,
        distance = dist_mode,
        similarity = sim_mode
    )
    return(res)
}

#' Compute spatial gradient field for input to DMT
#'
#' First, spatial gradients for each embedding dimension are computed for every
#' points (cell) by looking at each cell's neighbors. These point gradients can
#' be smoothed using bilateral/anisotropic filtering. Then gradients for triangles
#' are computed as the average of their vertices (3 for proper triangles, 2 for
#' degenerate exterior triangles at the boundaries). Finally, primal and dual
#' edge gradients are computed as the sum of the gradients for the two points and
#' two triangles, respectively, that are associated with each edge.
#'
#' @param dmt A list containing the mesh data structures:
#' * `pts` is a `N` x `2+M` data table with columns `X` and `Y` containing the
#'   coordinates of cells and additional metadata.
#' * `tris` is a `F` x `4` data table containing the X,Y coordinates of each
#'    triangle's centroid in the first two columns, and area and
#'    largest height of each triangle in the last two columns.
#' * `edges` is a `E` x `14` data table with columns `from_pt`, `to_pt`, `from_tri`, `to_tri`,
#'   `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`, `x1_tri`, `y0_tri`, `y1_tri`,
#'   `length_pt`, `length_tri`. If only one triangle uses an edge, then the `from_tri`,
#'   `x0_tri`, and `y0_tri` fields will contain NaN values.
#' * `tri_to_pt` is a `F` x `N` sparse matrix with value 1 at (i,j) if
#'   triangle i uses point j as a vertex.
#' * `udv_cells` contains cell embeddings stored in `embeddings`
#'   (a `N` x `D` matrix with `D`-dimensional embeddings for each cell)
#'   and `loadings` (a `G` x `D` matrix with gene loadings for each latent variable).
#' @param smooth_distance One of `c('none', 'euclidean', 'projected', 'constant')`.
#'   If either `smooth_distance` or `smooth_similarity` is `'none'` (the default),
#'   then no smoothing of the gradient field is conducted.
#' @param smooth_similarity One of `c('none', 'euclidean', 'projected', 'constant')`.
#'   If either `smooth_distance` or `smooth_similarity` is `'none'` (the default),
#'   then no smoothing of the gradient field is conducted.
#' @param smooth_iter Number of rounds of gradient smoothing.
#'
#' @returns A gradient field with the following attributes:
#' \item{pts}{A `2` x `D` x `N` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` latent variables at every point in space.}
#' \item{tris}{A `2` x `D` x `F` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` latent variables at every triangle in the mesh.
#'   Average of the vertices (3 for full triangles, 2 for degenerate triangles).}
#' \item{edges_pts}{A `2` x `D` x `E` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` latent variables at every primal edge (point-to-point) in the mesh.
#'   Average of the two endpoints.}
#' \item{edges_tris}{A `2` x `D` x `E` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` latent variables at every dual edge (triangle-to-triangle) in the mesh.
#'   Average of the two adjacent triangles.}
#'
#' @export
compute_gradients = function(
    dmt,
    smooth_distance='none', smooth_similarity='none', smooth_iter=1,
    on_edges = FALSE
) {
    if (on_edges) {
        return(compute_gradients_edges(
            dmt,
            smooth_distance=smooth_distance,
            smooth_similarity=smooth_similarity,
            smooth_iter=smooth_iter
        ))
    }

    field = list()

    ## First on points
    coords = as.matrix(dmt$pts[, .(X, Y)])
    adj = as.matrix(igraph::from_edgelist()$fun(as.matrix(dmt$edges[, .(from_pt, to_pt)]), directed=FALSE))
    # adj = as.matrix(igraph::from_edgelist()$fun(as.matrix(dmt$edges[!is.na(from_pt), .(from_pt, to_pt)]), directed=FALSE))
    embeddings = dmt$udv_cells$embeddings
    field$pts = estimate_field(coords, adj, embeddings)

    ## Smooth pt gradient first
    if (smooth_distance != 'none' & smooth_similarity != 'none') {
        for (i in seq_len(smooth_iter)) {
            field$pts = smooth_field(coords, field=field$pts, adj, include_self=TRUE, distance=smooth_distance, similarity=smooth_similarity)
        }
    }

    ## Field on tris
    ## Now, we need to do some tensor math to derive fields for each triangle.
    field$tris = array(dim = c(2, dim(field$pts)[2], nrow(dmt$tris)))
    field$tris[1, , ] = as.matrix(field$pts[1, , ] %*% Matrix::t(dmt$tri_to_pt))
    field$tris[2, , ] = as.matrix(field$pts[2, , ] %*% Matrix::t(dmt$tri_to_pt))
    field$tris = sweep(field$tris, 3, Matrix::rowSums(dmt$tri_to_pt), '/')

    ## Finally on edges
    # i_edges_prim = which(!is.na(dmt$edges$from_pt))
    # field_edges_prim = array(dim = c(2, dim(field_pts)[2], length(i_edges_prim)))
    field$edges_pts = array(dim = c(2, dim(field$pts)[2], nrow(dmt$edges)))
    y1 = factor(dmt$edges$from_pt, levels=1:nrow(dmt$pts))
    y2 = factor(dmt$edges$to_pt, levels=1:nrow(dmt$pts))
    adj_edges_to_pts = Matrix::sparse.model.matrix(~0+y1) + Matrix::sparse.model.matrix(~0+y2)
    field$edges_pts[1, , ] = as.matrix(field$pts[1, , ] %*% Matrix::t(adj_edges_to_pts))
    field$edges_pts[2, , ] = as.matrix(field$pts[2, , ] %*% Matrix::t(adj_edges_to_pts))
    field$edges_pts = field$edges_pts / 2   # divide by 2 to keep same units as points and triangles

    field$edges_tris = array(dim = c(2, dim(field$pts)[2], nrow(dmt$edges)))
    y1 = factor(dmt$edges$from_tri, levels=1:nrow(dmt$tris))
    y2 = factor(dmt$edges$to_tri, levels=1:nrow(dmt$tris))
    adj_edges_to_tris = Matrix::sparse.model.matrix(~0+y1) + Matrix::sparse.model.matrix(~0+y2)
    field$edges_tris[1, , ] = as.matrix(field$tris[1, , ] %*% Matrix::t(adj_edges_to_tris))
    field$edges_tris[2, , ] = as.matrix(field$tris[2, , ] %*% Matrix::t(adj_edges_to_tris))
    field$edges_tris = field$edges_tris / 2   # divide by 2 to keep same units as points and triangles

    return(field)
}

#' Compress a gradient field using SVD
#'
#' Expresses the `2` x `D` total derivative at each location as
#' a pair of `2`-dimensional vectors in the gradient and orthogonal
#' directions.
#'
#' @param field A gradient field with the following attributes:
#'   * `pts`: A `2` x `D` x `N` array in column-major ordering
#'     containing the spatial gradient in expression for each of
#'     `D` latent variables at every point in space.
#'   * `tris`: A `2` x `D` x `F` array in column-major ordering
#'     containing the spatial gradient in expression for each of
#'     `D` latent variables at every triangle in the mesh.
#'     Average of the vertices (3 for full triangles, 2 for degenerate triangles).
#'   * `edges_pts`: A `2` x `D` x `E` array in column-major ordering
#'     containing the spatial gradient in expression for each of
#'     `D` latent variables at every primal edge (point-to-point) in the mesh.
#'     Sum of the two endpoints.
#'   * `edges_tris`: A `2` x `D` x `E` array in column-major ordering
#'     containing the spatial gradient in expression for each of
#'     `D` latent variables at every dual edge (triangle-to-triangle) in the mesh.
#'     Sum of the two adjacent triangles.
#'
#' @returns A gradient field with the same attributes as the input, as well as compressed
#'   representations `pts_svd`, `tris_svd`, `edges_pts_svd`, and `edges_tris_svd`. Each of these
#'   is a `N` x `6` matrix with the following columns for
#'   each location:
#'   \item{dx_grad,dy_grad}{x,y directions of unit vector in the
#'     direction of greatest change (first singular vector).}
#'   \item{dx_ortho,dy_ortho}{x,y directions of unit vector orthogonal
#'     to the direction of greatest change (second singular vector).}
#'   \item{len_grad,len_ortho}{Magnitude of directional derivative in the
#'     gradient and orthogonal directions (singular values).}
#'
#' @export
compress_gradients_svd = function(field) {
    ## compress gradients into 2D representations
    field$pts_svd = compress_field_cpp(field$pts)
    field$tris_svd = compress_field_cpp(field$tris)
    field$edges_pts_svd = compress_field_cpp(field$edges_pts)
    field$edges_tris_svd = compress_field_cpp(field$edges_tris)

    colnames(field$pts_svd) = c('dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho', 'len_grad', 'len_ortho')
    colnames(field$tris_svd) = c('dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho', 'len_grad', 'len_ortho')
    colnames(field$edges_pts_svd) = c('dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho', 'len_grad', 'len_ortho')
    colnames(field$edges_tris_svd) = c('dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho', 'len_grad', 'len_ortho')

    return(field)
}

#' Compute a spatial gradient field at each point (cell)
#'
#' Distance between neighboring cells is normalized to unit distance
#' so that only the direction from each cell to its neighbors matters.
#' The gradient is then the average gradient in expression of each
#' embedding dimension between the index cell and its neighbors.
#'
#' @param coords A `N` x `2` matrix of cell coordinates.
#' @param embeddings A `N` x `D` matrix of cell embeddings.
#' @param adj A `N` x `N` sparse adjacency matrix
#'   in dgCMatrix format.
#'
#' @returns A `2` x `D` x `N` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` embedding dimensions at every point in space.
#'
#' @export
estimate_field = function(coords, adj, embeddings) {
    diag(adj) = 0 ## cannot include self by definition
    res = estimate_field_cpp(
        coords = coords,
        embeddings = embeddings,
        adj_i = adj@i,
        adj_p = adj@p
    )
    return(res)
}

# estimate_field = function(coords, adj, embeddings, include_self=TRUE) {
#     if (include_self) diag(adj) = 1
#     res = estimate_field_cpp(
#         coords = coords,
#         embeddings = embeddings,
#         adj_i = adj@i,
#         adj_p = adj@p
#     )

#     colnames(res) = c('len_grad', 'len_ortho', 'dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho')
#     res = data.table(res)
#     return(res)
# }

# smooth_field_old = function(coords, field, adj, embeddings, include_self=TRUE, reorient_grad_pos_y=TRUE) {
#     if (include_self) diag(adj) = 1
#     res = smooth_field_cpp(
#         pvec = diff(adj@p),
#         jvec = adj@i, ## keep it 0-indexed
#         field_coords = as.matrix(field[, c('dx_grad', 'dy_grad')] * field$len_grad),
#         coords = coords,
#         embeddings = embeddings,
#         adj_p = adj@p,
#         reorient_grad_pos_y
#     )

#     colnames(res) = c('len_grad', 'len_ortho', 'dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho')
#     res = data.table(res)
#     return(res)
# }

