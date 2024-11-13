#' Get separatrices that separate points into components with strong boundaries
#' 
#' Finds the the collection of edges that lie along separatrices. These are found by
#' tracing paths from saddles to critical triangles (maxima) in the dual spanning forest.
#' Because the dual spanning forest is a maximum spanning forest on triangles in the mesh,
#' these paths will follow ridges, separating points into components with the strongest boundaries.
#' 
#' Saddle edges are found by interesting the saddle edges found in the primal forest with
#' saddles edges found in the dual forest. All boundary edges in the dual forest are also
#' included as possible saddle edges.
#'
#' @param dmt A DMT data structure with the following attributes:
#'  * `edges`: Data structure of edges in the mesh.
#'  * `prim`: Data structure of primal minimum spanning forest.
#'  * `dual`: Data structure of dual maximum spanning forest.
#'
#' @returns A length `num_sep_edges` vector of edges (1-indexed) that make up the separatrices
#' that separate points into different components.
#' @export
dmt_get_separatrices = function(dmt) {
    ## Reduce graphs with intersecting saddles
    dual_saddles = c(
        dmt$dual$saddles, 
        which(dmt$edges$boundary==TRUE) ## let all boundaries be dual saddles
    )
    saddles = intersect(dmt$prim$saddles, dual_saddles)
    epaths = trace_paths(dmt, dmt$dual, saddles)        
    e_sep = get_e_sep(epaths, saddles, nrow(dmt$edges))    
    e_sep = as.integer(e_sep + 1)
    return(e_sep)
}

#' Set DMT scalar field values as the Frobenius norm of the total derivative
#' 
#' @inheritParams compute_gradients
#' @param field A gradient field with the compressed representations `pts_svd`,
#'   `tris_svd`, `edges_pts_svd`, and `edges_tris_svd`. Each of these
#'   is a `N` x `6` matrix with the following columns for each location:
#'   * `dx_grad,dy_grad`: x,y directions of unit vector in the
#'     direction of greatest change (first singular vector).
#'   * `dx_ortho,dy_ortho`: x,y directions of unit vector orthogonal
#'     to the direction of greatest change (second singular vector).
#'   * `len_grad,len_ortho`: Magnitude of directional derivative in the
#'     gradient and orthogonal directions (singular values).
#' 
#' @returns `dmt` with the following additional attributes:
#'   `dmt$pts$f`, `dmt$tris$f`, `dmt$edges$f_prim`, and `dmt$edges$f_dual`
#' 
#' @export
dmt_set_f = function(dmt, field) {
    dmt$pts$f = rowSums(field$pts_svd[, 5:6])
    dmt$tris$f = rowSums(field$tris_svd[, 5:6])
    dmt$edges$f_prim = rowSums(field$edges_pts_svd[, 5:6])
    dmt$edges$f_dual = rowSums(field$edges_tris_svd[, 5:6])
    return(dmt)
}

#' Construct primal minimum spanning forest on points
#'
#' Constructs a directed minimum spanning forest from point and edge scalar values
#' using a version of Prim's algorithm. Critical points (local minimum; or, more precisely,
#' an endpoint of an edge that is a local minimum) are used as roots for each tree in the
#' forest, and edges that would bridge two trees with different critical point roots are
#' marked as possible saddle edges.
#'
#' @param dmt A DMT data structure with the following attributes:
#'  * `pts$f`: Scalar field values defined for each point.
#'  * `edges$from_pt`: Index of first point joined by each edge.
#'  * `edges$to_pt`: Index of second point joined by each edge.
#'  * `edges$f_prim`: Scalar field values defined for each edge that connects two points.
#'
#' @returns A List with the following attributes (all indices are 1-indexed):
#'  * `edges`: A data.table with `forest_size` rows, where each row is a directed edge
#'    in the minimum spanning forest. There are six columns:
#'      * `from,to`: Index of source and target points for each edge.
#'      * `x0,y0`: Coordinates of source point for each edge.
#'      * `x1,y1`: Coordinates of target point for each edge.
#'  * `saddles`: A length `num_saddles` vector with edge indices for possible saddle edges.
#'  * `labels`: A length `num_points` vector of labels for the connected components in the
#'    minimum spanning tree. Each connected component is labeled by the index of its critical point.
#'  * `minima`: A length `num_critpts` vector of critical points (minima).
#'  * `parent`: A length `num_points` vector containing the parent (source) point for each
#'    point in the directed spanning forest. Critical points have no parent, so the value should be ignored.
#'  * `parent_edge`: A length `num_points` vector containing the directed edge that has
#'    each point as a target node. Critical points have no parent edge, so the value should be ignored.
#'
#' @export
do_primary_forest = function(dmt) {
    res = do_dmt_forest_cpp(
        -dmt$pts$f, 
        dmt$edges$from_pt - 1,
        dmt$edges$to_pt - 1,
        -dmt$edges$f_prim
    )
    names(res) = c('edges', 'saddles', 'labels', 'minima', 'parent', 'parent_edge')
    res$edges = data.table(res$edges)
    colnames(res$edges) = c('from', 'to')
    res$labels = as.integer(res$labels)    
    res$edges[, `:=`(x0 = dmt$pts$X[from], y0 = dmt$pts$Y[from], x1 = dmt$pts$X[to], y1 = dmt$pts$Y[to])]
    return(res)
}

#' Construct dual maximum spanning forest on triangles
#'
#' Constructs a directed maximum spanning forest from triangle and edge scalar values
#' using a version of Prim's algorithm. Critical triangles (local maximum; or, more precisely,
#' a triangle that contains an edge that is a local maximum) are used as roots for each tree in the
#' forest, and edges that would bridge two trees with different critical triangle roots are
#' marked as possible saddle edges.
#'
#' @param dmt A DMT data structure with the following attributes:
#'  * `tris$f`: Scalar field values defined for each triangle.
#'  * `edges$from_tri`: Index of first triangle joined by each edge.
#'  * `edges$to_tri`: Index of second triangle joined by each edge.
#'  * `edges$f_dual`: Scalar field values defined for each edge that connects two triangles.
#'
#' @returns A List with the following attributes (all indices are 1-indexed):
#'  * `edges`: A data.table with `forest_size` rows, where each row is a directed edge
#'    in the maximum spanning forest. There are six columns:
#'      * `from,to`: Index of source and target triangles for each edge.
#'      * `x0,y0`: Coordinates of source triangle for each edge.
#'      * `x1,y1`: Coordinates of target triangle for each edge.
#'  * `saddles`: A length `num_saddles` vector with edge indices for possible saddle edges.
#'  * `labels`: A length `num_triangles` vector of labels for the connected components in the
#'    maximum spanning tree. Each connected component is labeled by the index of its critical triangle.
#'  * `maxima`: A length `num_critpts` vector of critical triangles (maxima).
#'  * `parent`: A length `num_triangles` vector containing the parent (source) triangle for each
#'    triangle in the directed spanning forest. Critical triangles have no parent, so the value should be ignored.
#'  * `parent_edge`: A length `num_triangles` vector containing the directed edge that has
#'    each triangle as a target node. Critical triangles have no parent edge, so the value should be ignored.
#'
#' @export
do_dual_forest = function(dmt) {
    res = do_dmt_forest_cpp(
        dmt$tris$f, 
        dmt$edges$from_tri - 1,
        dmt$edges$to_tri - 1,
        dmt$edges$f_dual
    )
    names(res) = c('edges', 'saddles', 'labels', 'maxima', 'parent', 'parent_edge')
    res$edges = data.table(res$edges)
    colnames(res$edges) = c('from', 'to')
    res$labels = as.integer(res$labels)    
    res$edges[, `:=`(x0 = dmt$tris$X[from], y0 = dmt$tris$Y[from], x1 = dmt$tris$X[to], y1 = dmt$tris$Y[to])]
    return(res)
}

#' Trace all paths from saddles to dual critical points in the spanning forest.
#' 
#' Because the dual spanning forest is a maximum spanning forest on triangles in the mesh,
#' the paths that are traced between saddles and critical triangles (maxima) will follow
#' ridges, separating points into components with the strongest boundaries.
#' 
#' @param dmt A DMT data structure with the following attributes:
#'  * `edges$from_tri`: Index of first triangle joined by each edge (1-indexed).
#'  * `edges$to_tri`: Index of second triangle joined by each edge (1-indexed).
#' @param dual A dual maximum spanning forest on triangles:
#'  * `labels`: A length `num_triangles` vector of labels for the connected components in the
#'    maximum spanning tree. Each connected component is labeled by the index of its critical triangle.
#'  * `parent`: A length `num_triangles` vector containing the parent (source) triangle for each
#'    triangle in the directed spanning forest. Critical triangles have no parent, so the value should be ignored.
#'  * `parent_edge`: A length `num_triangles` vector containing the directed edge that has
#'    each triangle as a target node. Critical triangles have no parent edge, so the value should be ignored.
#' @param saddles A length `num_saddles` vector with edge indices for saddle edges (1-indexed).
#'
#' @returns A list of length `2*num_saddles` containing the two paths from each saddle edge to the
#'   two critical points that it joins. Each path is a numeric vector of edge indices (1-indexed).
#' @export
trace_paths = function(dmt, dual, saddles) {
    # saddles = intersect(prim$saddles, dual$saddles)    
    epaths = trace_epaths_cpp(
        saddles - 1, 
        dual$labels - 1, 
        dmt$edges$from_tri - 1, 
        dmt$edges$to_tri - 1, 
        dual$parent_edge - 1, 
        dual$parent - 1
    )
    return(epaths)
}


#' @export
get_boundary_shape = function(.SD, boundary) {
    ## Find endpoints of open line 
    i = as.integer(names(which(table(c(.SD$from_id, .SD$to_id)) == 1)))
    lvls = union(.SD$from_id, .SD$to_id) 
    .SD$from = as.integer(factor(.SD$from_id, lvls))
    .SD$to = as.integer(factor(.SD$to_id, lvls))
    

    ## First, find the incomplete path, from pt1 to pt2 
    pt1_i = c(.SD$from, .SD$to)[which(c(.SD$from_id, .SD$to_id) == i[1])]
    pt2_i = c(.SD$from, .SD$to)[which(c(.SD$from_id, .SD$to_id) == i[2])]
    path = igraph::from_edgelist()$fun(as.matrix(dplyr::select(.SD, from, to)), directed=FALSE) 
    path = igraph::shortest_paths(path,pt1_i, pt2_i)$vpath
    path = as.integer(path[[1]])
    path = dplyr::arrange(unique(bind_rows(
        dplyr::select(.SD, id = from, x=x0, y=y0),
        dplyr::select(.SD, id = to, x=x1, y=y1)
    )), id)[path, 2:3]
    

    ## All the steps below help to complete the cycle for this path 
    ## Find two closest points from incomplete shape to the boundary line 
    pt1 = sf::st_point(as.numeric(sf::st_drop_geometry(dmt$tris)[i[1], c('X', 'Y')]))
    pt2 = sf::st_point(as.numeric(sf::st_drop_geometry(dmt$tris)[i[2], c('X', 'Y')]))
    

    ## Split boundary into two lines using pt1 and pt2 
    lines = sf::st_split_line(boundary, pt1, pt2) ## this is our custom function 


    stopifnot(length(lines) == 2) ## must be true for closed linestrings 

    ## Try to close polygon with both lines
    ## Keep smaller area shape
    shapes = purrr::map(lines, function(line) {
        line = sf::st_coordinates(line)[, 1:2]
        ## order the line so that it starts at pt2 
        if (
            sum((head(line, 1) - sf::st_coordinates(pt2))^2) > 
            sum((tail(line, 1) - sf::st_coordinates(pt2))^2)
        ) {
            line = line[rev(seq_len(nrow(line))), ] ## reserve order 
        }
    
        colnames(line) = colnames(path)
        path = Reduce(rbind, list(path, line, path[1, ])) ## connect boundary line to dual forest path 
        shape = sf::st_polygon(list(as.matrix(path)))
        return(shape)
    })


    if (sf::st_area(shapes[[1]]) < sf::st_area(shapes[[2]])) {
        return(shapes[[1]])
    } else {
        return(shapes[[2]])
    }
}
    
#' Initialize tiles with shapes and other properties
#' 
#' @param dmt A DMT data structure with `pts`, `edges`, and `udv_cells` attributes.
#'   Assignment of points to tiles is provided in `dmt$pts$agg_id`, `dmt$edges$agg_from`,
#'   and `dmt$edges$agg_to`.
#' 
#' @returns A list data structure for tiles and their adjacencies. Includes the following attributes:
#'   * `meta_data`: A data table with attributes for tiles:
#'     * `X,Y`: Centroid of each tile.
#'     * `npts`: Number of points in each tile.
#'     * `shape`: A `sfc` list-column with the geometries for each tile.
#'     * `area,perimeter`: Area and perimeter of each tile.
#'   * `edges`: A data table with attributes for edges between adjacent tiles:
#'     * `from,to`: Tile IDs for the two tiles bordering this edge.
#'     * `x0,y0,x1,y1`: Centroid coordinates for the two tiles bordering this edge.
#'     * `area,npts`: Sum of areas and numbers of points in the two tiles bordering this edge.
#'     * `edge_length`: Total length of the border between the `from` and `to` tiles.
#'   * `pcs`: A `num_tiles` x `npcs` matrix with the average embedding value over all cells
#'     in each tile.
#' 
#' @export
dmt_init_tiles = function(dmt) {
    aggs = list()

    aggs$meta_data = dmt$pts[, .(X = mean(X), Y = mean(Y)), by = agg_id]
    aggs$meta_data = dplyr::rename(aggs$meta_data, id = agg_id)
    aggs$meta_data$npts = dmt$pts[, .N, by = agg_id][, N] 
    aggs$edges = unique(dmt$edges[, .(from=agg_from, to=agg_to)][
        from != to
    ][
        , `:=`(from = pmin(from, to), to = pmax(from, to))
    ])[
        , `:=`(
            x0 = aggs$meta_data$X[from],
            y0 = aggs$meta_data$Y[from],
            x1 = aggs$meta_data$X[to],
            y1 = aggs$meta_data$Y[to]
        )
    ][]
    
    aggs$meta_data$shape = trace_polygons(dmt, aggs)  
    aggs$meta_data[
        , `:=`(
            area = sf::st_area(shape), 
            perimeter = sf::st_length(sf::st_boundary(shape))
        )
    ]
    
    aggs$edges$area = aggs$meta_data$area[aggs$edges$from] + aggs$meta_data$area[aggs$edges$to]
    aggs$edges$npts = aggs$meta_data$npts[aggs$edges$from] + aggs$meta_data$npts[aggs$edges$to]
    
    ## compute edge length of edges 
    edge_lens = dmt$edges[
        agg_from != agg_to
    ][
        , `:=`(from = pmin(agg_from, agg_to), to = pmax(agg_from, agg_to))
    ][
        , len := sqrt((x0_tri-x1_tri)^2 + (y0_tri-y1_tri)^2)
    ][
        , .(edge_length = sum(len)) , by = .(from, to=to)
    ]
    
    aggs$edges = aggs$edges[edge_lens, on = c('from', 'to')]
    
    design = Matrix::sparse.model.matrix(~0+factor(agg_id), dmt$pts)
    design = design %*% Matrix::Diagonal(x = 1/Matrix::colSums(design)) 
    aggs$pcs = as.matrix(Matrix::t(design) %*% as.matrix(dmt$udv_cells$embeddings))
    return(aggs)
}

#' Assign points to tiles after DMT
#' 
#' Tiles are defined as points that are in the same connected component
#' in the mesh after discounting separatrix edges.
#' 
#' @param dmt A DMT data structure with `edges`, `pts`, and `e_sep` attributes.
#' 
#' @returns A DMT data structure with the following additional attributes:
#'   * `pts$agg_id`: Unique ID for the tile that each point belongs to.
#'   * `edges$agg_from`: Unique ID for the tile that `edges$from_pt` belongs to.
#'   * `edges$agg_to`: Unique ID for the tile that `edges$to_pt` belongs to.
#' @export
dmt_assign_tiles = function(dmt) {
    e = setdiff(seq_len(nrow(dmt$edges)), dmt$e_sep)
    # igraph::components(igraph::from_edgelist()$fun(as.matrix(dmt$edges)[e, c("from_pt", "to_pt")], FALSE))$membership
    g = Matrix::sparseMatrix(i = dmt$edges$from_pt[e], j = dmt$edges$to_pt[e], x = 1, dims = c(nrow(dmt$pts), nrow(dmt$pts)))
    g = igraph::from_adjacency()$fun(g, 'undirected')
    dmt$pts$agg_id = igraph::components(g)$membership
    dmt$edges$agg_from = dmt$pts$agg_id[dmt$edges$from_pt]
    dmt$edges$agg_to = dmt$pts$agg_id[dmt$edges$to_pt]    
    return(dmt)
}
