#' @export
GetTiles = function(
    X, Y, counts, 
    
    meta_data = NULL, meta_vars_include = NULL, 

    ###### STEP 0 ######
    npcs = 10, 
    ## Graph pruning
    prune_thresh_quantile = 0.95, 
    prune_min_cells = 10, 

    ###### STEP 1: GRADIENTS ######
    smooth_distance = c('none', 'euclidean', 'projected', 'constant')[1], 
    smooth_similarity = c('none', 'euclidean', 'projected', 'constant')[1], 

    ###### STEP 2: DMT ######

    ###### STEP 3: AGGREGATION ######
    max_npts = 50, 
    min_npts = 5, 
    alpha = 1, ## 0.2 = conservative merging, 2 = liberal merging 

    verbose = TRUE

) {

    ## STEP 0: PREPARE DATA STRUCTURES
    if (verbose) message('STEP 0: PREPARE DATA STRUCTURES')
    dmt = init_data(X, Y, counts, meta_data, meta_vars_include)
    dmt = prune_graph(dmt, thresh_quantile = prune_thresh_quantile, mincells = prune_min_cells) 
    dmt = add_exterior_triangles(dmt)

    dmt$udv_cells = do_pca(dmt$counts, npcs)

    ## STEP 1: GRADIENTS 
    if (verbose) message('STEP 1: GRADIENTS ')
    field = compute_gradients(dmt, smooth_distance, smooth_similarity)
    field = compress_gradients_svd(field)    
    
    ## STEP 2: DMT
    if (verbose) message('STEP 2: DMT')
    dmt = dmt_set_f(dmt, field)    
    dmt$prim = do_primary_forest(dmt)
    dmt$dual = do_dual_forest(dmt)
    dmt$e_sep = dmt_get_separatrices(dmt)    
    dmt = dmt_assign_tiles(dmt)
    aggs = dmt_init_tiles(dmt)    
    
    ## STEP 3: AGGREGATION
    if (verbose) message('STEP 3: AGGREGATION')
    ## First, main aggregation
    aggs = init_scores(aggs, agg_mode=2, alpha=alpha, max_npts=max_npts)
    aggs = merge_aggs(aggs, agg_mode=2, max_npts=max_npts) 
    dmt = update_dmt_aggid(dmt, aggs)
    aggs = update_agg_shapes(dmt, aggs)    

    ## Then, clean up stray small aggs
    aggs = init_scores(aggs, agg_mode=3, alpha=alpha, min_npts=min_npts)  
    aggs = merge_aggs(aggs, agg_mode=3, min_npts=min_npts)
    dmt = update_dmt_aggid(dmt, aggs)
    aggs = update_agg_shapes(dmt, aggs)    
    
    return(list(dmt=dmt, aggs=aggs))
}


