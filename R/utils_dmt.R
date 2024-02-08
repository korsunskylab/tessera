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

#' @export
dmt_set_f = function(dmt, field) {
    dmt$pts$f = rowSums(field$pts_svd[, 5:6])
    dmt$tris$f = rowSums(field$tris_svd[, 5:6])
    dmt$edges$f_prim = rowSums(field$edges_pts_svd[, 5:6])
    dmt$edges$f_dual = rowSums(field$edges_tris_svd[, 5:6])
    return(dmt)
}


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
    path = igraph::from_edgelist()$fun(as.matrix(dplyr::select(.SD, from, to)), directed=FALSE) %>% 
        igraph::shortest_paths(pt1_i, pt2_i) %>% with(vpath)
    path = as.integer(path[[1]])
    path = dplyr::arrange(unique(bind_rows(
        dplyr::select(.SD, id = from, x=x0, y=y0),
        dplyr::select(.SD, id = to, x=x1, y=y1)
    )), id)[path, 2:3]
    

    ## All the steps below help to complete the cycle for this path 
    ## Find two closest points from incomplete shape to the boundary line 
    pt1 = st_point(as.numeric(st_drop_geometry(dmt$tris)[i[1], c('X', 'Y')]))
    pt2 = st_point(as.numeric(st_drop_geometry(dmt$tris)[i[2], c('X', 'Y')]))
    

    ## Split boundary into two lines using pt1 and pt2 
    lines = st_split_line(boundary, pt1, pt2) ## this is our custom function 


    stopifnot(length(lines) == 2) ## must be true for closed linestrings 

    ## Try to close polygon with both lines
    ## Keep smaller area shape
    shapes = purrr::map(lines, function(line) {
        line = st_coordinates(line)[, 1:2]
        ## order the line so that it starts at pt2 
        if (
            sum((head(line, 1) - st_coordinates(pt2))^2) > 
            sum((tail(line, 1) - st_coordinates(pt2))^2)
        ) {
            line = line[rev(seq_len(nrow(line))), ] ## reserve order 
        }
    
        colnames(line) = colnames(path)
        path = Reduce(rbind, list(path, line, path[1, ])) ## connect boundary line to dual forest path 
        shape = st_polygon(list(as.matrix(path)))
        return(shape)
    })


    if (st_area(shapes[[1]]) < st_area(shapes[[2]])) {
        return(shapes[[1]])
    } else {
        return(shapes[[2]])
    }
}
    

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
            area = st_area(shape), 
            perimeter = st_length(st_boundary(shape))
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
