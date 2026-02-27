#' Run Tessera spatial segmentation pipeline
#'
#' S3 generic that dispatches to [RunTessera.default()] for matrix/vector
#' inputs and [RunTessera.Seurat()] for Seurat objects.
#'
#' @family RunTessera
#' @rdname RunTessera
#' @inheritDotParams RunTessera.default -X -Y -counts -embeddings -loadings -meta_data
#'
#' @returns See [RunTessera.default()].
#'
#' @export
RunTessera = function(...) {
	UseMethod("RunTessera")
}


#' Applies Tessera on a Seurat object
#'
#' @family RunTessera
#' @rdname RunTessera.Seurat
#'
#' @param obj Seurat object with spatial coordinates (and optionally,
#'   pre-computed cell embeddings) stored as dimensional reductions.
#' @param spatial Name of dimensional reduction where the cells' x/y
#'   coordinates are stored.
#' @param embeddings Name of dimensional reduction where pre-computed
#'   single-cell embeddings are stored. If NULL, PCA is computed from counts.
#'   If provided, `npcs` is ignored.
#' @param assay Seurat assay to use for counts. Defaults to `DefaultAssay`.
#' @param raw_results If TRUE, return the raw list from [RunTessera.default()].
#' @param tile.id.name Column name for tile IDs in the input Seurat object.
#' @param reduction.name Name of the tile-level PCA reduction in the tile Seurat object.
#' @param graph.name Name of the tile adjacency graph slot in the tile Seurat object.
#' @inheritDotParams RunTessera.default -X -Y -counts -embeddings -loadings -meta_data
#'
#' @returns A List with:
#' * `obj`: input Seurat object updated with tile assignments in `meta.data`.
#' * `tile_obj`: Seurat object with one item per tile.
#'
#' @export
RunTessera.Seurat = function(
	obj,
	spatial,
	embeddings    = NULL,
	assay         = NULL,
	raw_results   = FALSE,
	tile.id.name  = 'tile_id',
	reduction.name = 'pca',
	graph.name    = 'tile_adj',
	...
) {
	rlang::check_installed("Seurat", reason = "to use `RunTessera.Seurat()`")

	if ((!is.null(embeddings)) && (embeddings != 'harmony')) {
		warning('It is recommended to use harmonized cell embeddings as input for Tessera')
	}

	if (is.null(assay)) {
		assay = Seurat::DefaultAssay(obj)
	}

	if (!is.null(embeddings)) {
		emb  = Seurat::Embeddings(obj, reduction = embeddings)
		load = Seurat::Loadings(obj, reduction = embeddings)
	} else {
		emb  = NULL
		load = NULL
	}

	res = RunTessera(
		X          = Seurat::Embeddings(obj, spatial)[, 1],
		Y          = Seurat::Embeddings(obj, spatial)[, 2],
		counts     = obj[[assay]]$counts,
		embeddings = emb,
		loadings   = load,
		meta_data  = obj@meta.data,
		...
	)

	if (raw_results) return(res)

	stopifnot(all(colnames(res$aggs$counts) == res$aggs$meta_data$id))

	# Add tile IDs to input Seurat object
	if (tile.id.name %in% colnames(obj@meta.data)) {
		tile.id.name = tail(make.unique(c(colnames(obj@meta.data), tile.id.name)), n = 1)
		warning(paste0(
			'To avoid overwriting existing meta.data variables, tile.id.name is set to ',
			tile.id.name
		))
	}
	obj@meta.data[[tile.id.name]] = as.character(NA)
	obj@meta.data[[tile.id.name]][res$dmt$pts$ORIG_ID] = res$dmt$pts$agg_id
	obj@meta.data[[tile.id.name]] = factor(obj@meta.data[[tile.id.name]])

	# Build tile-level Seurat object
	meta.data = data.frame(dplyr::select(res$aggs$meta_data, -shape))
	row.names(meta.data) = meta.data$id
	stopifnot(all(colnames(res$aggs$counts) == res$aggs$meta_data$id))
	stopifnot(all(colnames(res$aggs$counts) == row.names(meta.data)))

	tile_obj = Seurat::CreateSeuratObject(
		counts    = res$aggs$counts,
		meta.data = meta.data
	)
	tile_obj@meta.data$shape = res$aggs$meta_data$shape
	row.names(res$aggs$pcs)  = colnames(tile_obj)
	tile_obj[[reduction.name]] = Seurat::CreateDimReducObject(
		embeddings = res$aggs$pcs, key = reduction.name
	)
	tile_obj[[graph.name]] = Seurat::as.Graph(res$aggs$adj)

	list(obj = obj, tile_obj = tile_obj)
}


#' Run full Tessera segmentation pipeline on spatial transcriptomics data
#'
#' Orchestrates four pipeline stages: (1) Delaunay mesh construction and
#' pruning, (2) spatial gradient field estimation with optional bilateral
#' smoothing, (3) Discrete Morse Theory segmentation, and (4) agglomerative
#' tile merging.
#'
#' Internally calls [make_cells()], [make_mesh()], [compute_field()],
#' [compute_morse()], [run_dmt()], and [make_tiles()].
#'
#' @family RunTessera
#' @rdname RunTessera.default
#'
#' @param X,Y Numeric vectors of length `num_cells` with cell coordinates.
#' @param counts A `num_genes` x `num_cells` gene-by-cell count matrix.
#'   Optional if `embeddings` are provided.
#' @param embeddings A `num_cells` x `num_dim` matrix of pre-computed cell
#'   embeddings. If NULL, PCA is computed from counts. If provided, `npcs`
#'   is ignored.
#' @param loadings Optional `num_genes` x `num_dim` matrix of gene loadings.
#' @param meta_data Data frame of per-cell metadata (one row per cell).
#' @param meta_vars_include Column names from `meta_data` to attach to each
#'   surviving cell in the output `dmt$pts` table.
#' @param group.by Column in `meta_data` identifying distinct samples or FOVs.
#'   Tessera is run separately for each group. If NULL, all cells are treated
#'   as a single group.
#' @param npcs Number of PCs to compute when `embeddings` is NULL.
#' @param prune_thresh_quantile Quantile of edge lengths above which edges are
#'   pruned. Default 0.95.
#' @param prune_min_cells Minimum cells per connected component after pruning.
#'   Default 10.
#' @param prune_thresh Absolute edge-length threshold; overrides
#'   `prune_thresh_quantile` when not NA. Default NA.
#' @param smooth_distance,smooth_similarity Distance and similarity weighting
#'   for bilateral gradient smoothing. Each is one of `'none'`, `'euclidean'`,
#'   `'projected'`, or `'constant'`. Default `'projected'`.
#' @param smooth_iter Number of bilateral smoothing iterations. Default 1.
#' @param max_npts Maximum cells per tile in the main agglomeration pass.
#'   Default 50.
#' @param min_npts Minimum cells per tile in the cleanup pass. Default 5.
#' @param alpha GMM sigma multiplier for transcriptional similarity scoring.
#'   `0.2` = conservative merging, `2` = liberal merging. Default 1.
#' @param .progress Show `furrr` progress bar. Default TRUE.
#' @param .options `furrr_options` object passed to [furrr::future_map()].
#'   Defaults to `furrr_options(stdout=TRUE, seed=TRUE)`.
#' @param future.globals.maxSize Maximum size in bytes of globals exported to
#'   parallel workers. Default 8 GB.
#' @param consolidate If TRUE (default), combines results across all groups
#'   into a single `list(dmt, aggs)`. If FALSE, returns a named list with one
#'   entry per group.
#' @param verbose Print per-stage progress messages. Default FALSE.
#'
#' @returns If `consolidate == TRUE` (default), a list with:
#' \item{dmt}{Cell-level mesh data structure with slots `pts`, `tris`,
#'   `edges`, `tri_to_pt`, `counts`, `udv_cells`, `prim`, `dual`, `e_sep`.
#'   `dmt$pts` has columns `X`, `Y`, `ORIG_ID`, `f`, `agg_id`, plus any
#'   columns in `meta_vars_include`.}
#' \item{aggs}{Tile-level data structure with slots `meta_data`, `counts`,
#'   `pcs`, `edges`, `adj`.}
#'
#' If `consolidate == FALSE`, a named list with one `list(dmt, aggs)` per
#' group.
#'
#' @export
RunTessera.default = function(
	X, Y,
	counts    = NULL,
	embeddings = NULL,
	loadings  = NULL,
	meta_data = NULL,
	meta_vars_include = NULL,
	group.by  = NULL,

	###### STEP 0 ######
	npcs = 10,
	## Graph pruning
	prune_thresh_quantile = 0.95,
	prune_min_cells       = 10,
	prune_thresh          = NA,

	###### STEP 1: GRADIENTS ######
	smooth_distance   = c('none', 'euclidean', 'projected', 'constant')[3],
	smooth_similarity = c('none', 'euclidean', 'projected', 'constant')[3],
	smooth_iter       = 1,

	###### STEP 2: DMT ######

	###### STEP 3: AGGREGATION ######
	max_npts = 50,
	min_npts = 5,
	alpha    = 1,

	## future_map parameters
	.progress = TRUE,
	.options  = NULL,
	future.globals.maxSize = 8 * 1024^3,  # 8 GB

	consolidate = TRUE,
	verbose     = FALSE
) {
	if (length(X) != length(Y))
		stop('X and Y have different lengths')

	if (is.null(counts)) {
		if (is.null(embeddings))
			stop('Both counts and embeddings are missing. Must supply at least one.')
		warning('No counts provided. Using placeholder values.')
		counts = Matrix::sparseMatrix(c(), c(), dims = c(0, length(X)))
	}
	if (is.null(embeddings))
		warning('No embeddings provided. Calculating embeddings using PCA.')

	if (is.null(group.by)) {
		warning('No value for group.by provided. Analyzing as a single sample.')
		if (is.null(meta_data)) {
			group.by  = 'group'
			meta_data = data.frame(group = factor(rep(1, length(X))))
		} else {
			group.by = tail(make.unique(c(colnames(meta_data), 'group')), n = 1)
			meta_data[[group.by]] = factor(rep(1, length(X)))
		}
	} else {
		if (!(group.by %in% colnames(meta_data))) {
			stop('group.by must be a column of meta_data')
		} else if (!(group.by %in% colnames(meta_vars_include))) {
			# NOTE: colnames(character_vector) == NULL, so this branch always
			# fires when meta_vars_include is a plain character vector.
			meta_vars_include = c(meta_vars_include, group.by)
		}
	}

	# Temporarily increase global size limit for future
	old_opts = options(future.globals.maxSize = future.globals.maxSize)
	on.exit(options(old_opts), add = TRUE)

	if (is.null(.options))
		.options = furrr::furrr_options(stdout = TRUE, seed = TRUE)

	groups        = unique(as.character(unique(meta_data[[group.by]])))
	names(groups) = groups

	res = furrr::future_map(groups, function(group) {

		idx = which(meta_data[[group.by]] == group)

		## ── STEP 0: BUILD MESH ────────────────────────────────────────────────
		if (verbose) message('STEP 0: PREPARE DATA STRUCTURES')

		group_counts = counts[, idx, drop = FALSE]
		group_meta   = data.table::as.data.table(meta_data)[idx, ]
		coords       = matrix(
			c(X[idx], Y[idx]), ncol = 2,
			dimnames = list(NULL, c('x', 'y'))
		)

		cells = make_cells(coords, NULL, group_counts, group_meta)
		mesh  = make_mesh(cells, list(
			prune_thresh_quantile = prune_thresh_quantile,
			prune_mincells        = prune_min_cells,
			prune_thresh          = prune_thresh
		))

		## ── PCA or user-supplied embeddings ───────────────────────────────────
		pruned_counts = group_counts[, mesh$pts$ORIG_ID, drop = FALSE]
		if (is.null(embeddings)) {
			udv          = do_pca(pruned_counts, npcs)
			grp_embed    = udv$embeddings
			grp_loadings = udv$loadings
		} else {
			grp_embed    = embeddings[idx, , drop = FALSE][
				as.integer(mesh$pts$ORIG_ID), , drop = FALSE
			]
			grp_loadings = loadings
		}
		cells$embeddings = list(pca = grp_embed)

		## ── STEP 1: GRADIENTS ─────────────────────────────────────────────────
		if (verbose) message('STEP 1: GRADIENTS')
		field = compute_field(cells, mesh, list(
			smooth_distance   = smooth_distance,
			smooth_similarity = smooth_similarity,
			smooth_iter       = smooth_iter
		))

		## ── STEP 2: DMT ───────────────────────────────────────────────────────
		if (verbose) message('STEP 2: DMT')
		mesh$morse = compute_morse(field, mesh)
		dmt_run    = run_dmt(mesh)

		## ── STEP 3: TILE EXTRACTION ──────────────────────────────────────────
		if (verbose) message('STEP 3: TILE EXTRACTION')
		tiles_init = make_tiles(cells, mesh, dmt_run)

		## ── STEP 4: TILE MERGING ─────────────────────────────────────────────
		if (verbose) message('STEP 4: TILE MERGING')
		tiles = merge_tiles(
			tiles_init,
			list(alpha = alpha, max_npts = max_npts, min_npts = min_npts)
		)

		## ── Backward-compatible dmt struct ────────────────────────────────────
		pts_compat = data.table::copy(mesh$pts)
		pts_compat[, f      := mesh$morse$f$pts]
		pts_compat[, agg_id := tiles$cell_ids$tile_id]

		# Attach meta_vars (includes group.by when a group.by was supplied)
		if (!is.null(meta_vars_include) && length(meta_vars_include) > 0) {
			for (v in meta_vars_include)
				pts_compat[[v]] = cells$meta[[v]][mesh$pts$ORIG_ID]
		}

		# Map ORIG_ID from within-group index to global cell index
		pts_compat[, ORIG_ID := idx[ORIG_ID]]

		edges_compat = data.table::copy(mesh$edges)
		edges_compat[, f_prim   := mesh$morse$f$edges_prim]
		edges_compat[, f_dual   := mesh$morse$f$edges_dual]
		edges_compat[, agg_from := pts_compat$agg_id[from_pt]]
		edges_compat[, agg_to   := pts_compat$agg_id[to_pt]]

		tris_compat = data.table::copy(mesh$triangles)
		tris_compat[, f := mesh$morse$f$tris]

		dmt_compat = list(
			pts       = pts_compat,
			tris      = tris_compat,
			edges     = edges_compat,
			tri_to_pt = mesh$tri_to_pt,
			counts    = pruned_counts,
			udv_cells = list(loadings = grp_loadings, embeddings = grp_embed),
			prim      = mesh$morse$prim,
			dual      = mesh$morse$dual,
			e_sep     = mesh$morse$e_sep
		)

		aggs_compat = list(
			meta_data = tiles$meta_data,
			counts    = tiles$counts,
			pcs       = tiles$pcs,
			edges     = tiles$edges,
			adj       = tiles$adj
		)

		# Set group labels
		stopifnot(all(dmt_compat$pts[[group.by]] == group))
		aggs_compat$meta_data[[group.by]] = group
		dmt_compat$pts[[group.by]]        = group
		aggs_compat$edges[[group.by]]     = group

		list(dmt = dmt_compat, aggs = aggs_compat)

	}, .progress = .progress, .options = .options)

	if (consolidate) {
		if (length(res) > 1)
			res = ConsolidateResults(res, group.by)
		else
			res = res[[1]]
		res$aggs = AddAggsAdjacencyMatrix(res$aggs)
		return(res)
	} else {
		return(res)
	}
}
