#' @export
smooth_field = function(coords, field, adj, include_self=TRUE,
                        distance="euclidean", similarity="euclidean") {
    if (include_self) diag(adj) = 1

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


#' @export
compute_gradients = function(dmt, smooth_distance=NULL, smooth_similarity=NULL) {
    field = list()

    ## First on points 
    coords = as.matrix(dmt$pts[, .(X, Y)])
    adj = as.matrix(igraph::from_edgelist()$fun(as.matrix(dmt$edges[, .(from_pt, to_pt)]), directed=FALSE))
    # adj = as.matrix(igraph::from_edgelist()$fun(as.matrix(dmt$edges[!is.na(from_pt), .(from_pt, to_pt)]), directed=FALSE))
    embeddings = dmt$udv_cells$embeddings
    field$pts = estimate_field(coords, adj, embeddings)

    ## Smooth pt gradient first
    # if (!is.null(smooth_distance) & !is.null(smooth_similarity))
    if (smooth_distance != 'none' & smooth_similarity != 'none')
        field$pts = smooth_field(coords, field=field$pts, adj, include_self=TRUE, distance=smooth_distance, similarity=smooth_similarity) 
    
    ## Field on tris
    ## Now, we need to do some tensor math to derive fields for each triangle. 
    field$tris = array(dim = c(2, dim(field$pts)[2], nrow(dmt$tris)))
    field$tris[1, , ] = as.matrix(field$pts[1, , ] %*% t(dmt$tri_to_pt))
    field$tris[2, , ] = as.matrix(field$pts[2, , ] %*% t(dmt$tri_to_pt))
    field$tris = sweep(field$tris, 3, rowSums(dmt$tri_to_pt), '/')

    ## Finally on edges 
    # i_edges_prim = which(!is.na(dmt$edges$from_pt))
    # field_edges_prim = array(dim = c(2, dim(field_pts)[2], length(i_edges_prim)))
    field$edges_pts = array(dim = c(2, dim(field$pts)[2], nrow(dmt$edges)))
    y1 = factor(dmt$edges$from_pt, levels=1:nrow(dmt$pts))
    y2 = factor(dmt$edges$to_pt, levels=1:nrow(dmt$pts))
    adj_edges_to_pts = Matrix::sparse.model.matrix(~0+y1) + Matrix::sparse.model.matrix(~0+y2)
    field$edges_pts[1, , ] = as.matrix(field$pts[1, , ] %*% t(adj_edges_to_pts))
    field$edges_pts[2, , ] = as.matrix(field$pts[2, , ] %*% t(adj_edges_to_pts))

    field$edges_tris = array(dim = c(2, dim(field$pts)[2], nrow(dmt$edges)))
    y1 = factor(dmt$edges$from_tri, levels=1:nrow(dmt$tris))
    y2 = factor(dmt$edges$to_tri, levels=1:nrow(dmt$tris))
    adj_edges_to_tris = Matrix::sparse.model.matrix(~0+y1) + Matrix::sparse.model.matrix(~0+y2)
    field$edges_tris[1, , ] = as.matrix(field$tris[1, , ] %*% t(adj_edges_to_tris))
    field$edges_tris[2, , ] = as.matrix(field$tris[2, , ] %*% t(adj_edges_to_tris))
    
    return(field)
}

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

