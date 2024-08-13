## All functions for performing tile-based analysis on multiple samples 

#' @export
dmt_bind2 = function(dmt1, dmt2) {
    return(list(
        pts = rbindlist(list(dmt1$pts, dmt2$pts)),
        edges = rbindlist(list(dmt1$edges, dmt2$edges))
    ))
}



#' @export
aggs_bind2 = function(aggs1, aggs2) {
    return(list(
        meta_data = rbindlist(list(aggs1$meta_data, aggs2$meta_data)), 
        edges = rbindlist(list(aggs1$edges, aggs2$edges)), 
        counts = Matrix::cbind2(aggs1$counts, aggs2$counts), 
        pcs = Matrix::rbind2(aggs1$pcs, aggs2$pcs)
    ))
}

#' @export
aggs_bind = function(objs) {
    ## Smarter order for binding should help with sparse matrix binding 
    aggs_list = purrr::map(objs, 'aggs')
    naggs = length(aggs_list)
    sizes = aggs_list %>% map('meta_data') %>% map_int(nrow)
    merge_next = order(sizes)[1:2]
    for (i in seq_len(naggs - 1)) {
        aggs_list[[merge_next[1]]] = aggs_bind2(aggs_list[[merge_next[1]]], aggs_list[[merge_next[2]]])
        aggs_list[[merge_next[2]]] = NULL
        sizes = aggs_list %>% map('meta_data') %>% map_int(nrow)
        merge_next = order(sizes)[1:2]
    }
    return(aggs_list[[1]])
}



#' @export
bind_objs = function(objs) {
    ## First, simplify and prepare objects for merging 
    objs = imap(objs, function(obj, sample_name) {
        aggs = obj$aggs
        dmt = obj$dmt
        
        ## Add sample_name as prefix to all IDs
        colnames(aggs$counts) = paste0(sample_name, '_', colnames(aggs$counts))
        aggs$meta_data$id = paste0(sample_name, '_', seq_len(nrow(aggs$meta_data)))
        aggs$meta_data$sample_id = sample_name
        aggs$edges$sample_id = sample_name 
        
        dmt$edges$sample_id = sample_name 
        dmt$pts$sample_id = sample_name
        dmt$pts$id = paste0(sample_name, '_', seq_len(nrow(dmt$pts)))
        dmt$pts$agg_id = paste0(sample_name, '_', dmt$pts$agg_id)

        return(list(
            aggs = aggs[c('meta_data', 'edges', 'counts')], 
            dmt = dmt[c('pts', 'edges')]
            # dmt = dmt[c('pts', 'counts')]
        ))
    })
    
    # aggs = objs %>% purrr::map('aggs') %>% purrr::reduce(aggs_bind2)
    aggs = aggs_bind(objs)

    ## Only save meta_data from points
    ## Assume that counts should be saved elsewhere
    dmt = objs %>% purrr::map('dmt') %>% purrr::reduce(dmt_bind2)# %>% with(pts)
    aggs$pts = dmt$pts
    aggs$pt_edges = dmt$edges
    
    ## Re-index edges between aggs
    aggs$edges$from = match(paste0(aggs$edges$sample_id, '_', aggs$edges$from), aggs$meta_data$id)
    aggs$edges$to = match(paste0(aggs$edges$sample_id, '_', aggs$edges$to), aggs$meta_data$id)

    ## Re-index edges between points 
    aggs$pt_edges$from_pt = match(paste0(aggs$pt_edges$sample_id, '_', aggs$pt_edges$from_pt), aggs$pts$id)
    aggs$pt_edges$to_pt = match(paste0(aggs$pt_edges$sample_id, '_', aggs$pt_edges$to_pt), aggs$pts$id)
    # aggs$pt_edges$from_tri = match(paste0(aggs$pt_edges$sample_id, '_', aggs$pt_edges$from_tri), aggs$pts$id)
    # aggs$pt_edges$to_tri = match(paste0(aggs$pt_edges$sample_id, '_', aggs$pt_edges$to_tri), aggs$pts$id)
    
    
    # return(aggs)
    return(list(aggs=aggs, dmt=dmt))
}
