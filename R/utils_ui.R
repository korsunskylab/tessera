


#' Consolidate Tessera results from multiple samples (groups) after constructing
#' Tessera tiles separately on cells from each group.
#' 
#' @param res Output of running RunTessera (when consolidate == FALSE).
#' @param group.by Name of metadata variable that identifies distinct groups.
#' 
#' @export
ConsolidateResults = function(res, group.by) {
    all_aggs = list()
    all_aggs$meta_data = rbindlist(lapply(names(res), function(group) {res[[group]]$aggs$meta_data}))
    all_aggs$meta_data[[group.by]] = factor(all_aggs$meta_data[[group.by]])
    all_aggs$meta_data$id = paste0(all_aggs$meta_data[[group.by]], '_', all_aggs$meta_data$id)
    if (length(unique(all_aggs$meta_data$id)) != length(all_aggs$meta_data$id)) {
        stop('not all tile IDs are unique')
    }

    all_aggs$counts = do.call(cbind, lapply(names(res), function(group) {res[[group]]$aggs$counts}))
    colnames(all_aggs$counts) = all_aggs$meta_data$id

    all_aggs$pcs = do.call(rbind, lapply(names(res), function(group) {res[[group]]$aggs$pcs}))

    all_aggs$edges = rbindlist(lapply(names(res), function(group) {res[[group]]$aggs$edges}))
    all_aggs$edges[[group.by]] = factor(all_aggs$edges[[group.by]])
    all_aggs$edges$from = paste0(all_aggs$edges[[group.by]], '_', all_aggs$edges$from)
    all_aggs$edges$to = paste0(all_aggs$edges[[group.by]], '_', all_aggs$edges$to)
    
    all_dmt = list()
    all_dmt$pts = rbindlist(lapply(names(res), function(group) {res[[group]]$dmt$pts}))
    all_dmt$pts[[group.by]] = factor(all_dmt$pts[[group.by]])
    all_dmt$pts$agg_id = paste0(all_dmt$pts[[group.by]], '_', all_dmt$pts$agg_id)

    return(list(dmt=all_dmt, aggs=all_aggs))
}

#' Construct tile adjacency matrix from consolidated RunTessera output.
#' 
#' @param aggs Aggregated tile information after consolidation.
#' 
#' @export
AddAggsAdjacencyMatrix = function(aggs) {
    aggs$adj <- as.matrix(igraph::graph_from_data_frame(
        d = data.frame(from = aggs$edges$from, to = aggs$edges$to),
        vertices = data.frame(name = aggs$meta_data$id),
        directed = FALSE))
    stopifnot(all(colnames(aggs$adj) == aggs$meta_data$id))
    stopifnot(all(rownames(aggs$adj) == aggs$meta_data$id))
    return(aggs)
}

