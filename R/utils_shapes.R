trace_polygons = function(dmt, aggs) {
    res = trace_polygons_cpp(
        as.matrix(dmt$edges), 
        nrow(aggs$meta_data),
        nrow(dmt$tris), 
        dmt$pts$agg_id
    ) 
    shapes = res %>% 
        purrr::map(list) %>% 
        purrr::map(st_polygon) %>% 
        st_sfc()
    shapes = st_make_valid(shapes)

    return(shapes)
}
