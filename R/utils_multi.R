## All functions for performing tile-based analysis on multiple samples 

#' @export
dmt_bind2 = function(dmt1, dmt2) {
    return(list(pts = bind_rows(dmt1$pts, dmt2$pts)))
    # return(list(
    #     pts = bind_rows(dmt1$pts, dmt2$pts)
    #     # counts = Matrix::cbind2(dmt1$counts, dmt2$counts)
    # ))
}

#' @export
aggs_bind2 = function(aggs1, aggs2) {
    return(list(
        meta_data = bind_rows(aggs1$meta_data, aggs2$meta_data), 
        edges = bind_rows(aggs1$edges, aggs2$edges), 
        counts = Matrix::cbind2(aggs1$counts, aggs2$counts)
    ))

}

#' @export
bind_objs = function(objs) {
    ## First, simplify and prepare objects for merging 
    objs = imap(objs, function(obj, sample_name) {
        aggs = obj$aggs
        dmt = obj$dmt

        colnames(aggs$counts) = paste0(sample_name, '_', colnames(aggs$counts))
        aggs$meta_data$id = paste0(sample_name, '_', aggs$meta_data$id)
        aggs$meta_data$sample_id = sample_name
        aggs$edges$sample_id = sample_name 
        dmt$pts$sample_id = sample_name
        dmt$pts$agg_id = paste0(sample_name, '_', dmt$pts$agg_id)

        return(list(
            aggs = aggs[c('meta_data', 'edges', 'counts')], 
            dmt = dmt[c('pts')]
            # dmt = dmt[c('pts', 'counts')]
        ))
    })
    
    
    # dmt = objs %>% purrr::map('dmt') %>% purrr::reduce(dmt_bind2)
    aggs = objs %>% purrr::map('aggs') %>% purrr::reduce(aggs_bind2)

    ## Only save meta_data from points
    ## Assume that counts should be saved elsewhere
    aggs$pts = objs %>% purrr::map('dmt') %>% purrr::reduce(dmt_bind2) %>% with(pts)
    return(aggs)
    # return(list(aggs=aggs, dmt=dmt))
}
