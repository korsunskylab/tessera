#' Contruct shapes that outline each tile
#' 
#' @param dmt A DMT mesh data structure with `pts`, `edges`, and `tris` attributes.
#'   Assignment of points to tiles is provided in `dmt$pts$agg_id`.
#' @param aggs A tile data structure used to specify the total number of tiles.
#' 
#' @returns A `sfc` list-column with the geometries for each tile.
#' 
#' @export
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
