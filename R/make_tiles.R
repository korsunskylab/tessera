#' Build tile adjacency matrix from a tiles list
#'
#' Constructs a symmetric sparse adjacency matrix from the edge list in
#' `tiles$edges`.  Tile IDs from `tiles$meta_data$id` are used as row/column
#' names.
#'
#' @param tiles A tiles list as returned by [make_tiles()].
#'
#' @returns A `num_tiles × num_tiles` symmetric dgCMatrix where entry (i, j)
#'   is 1 if tiles i and j share a border.
#'
#' @seealso [make_tiles()]
#' @export
make_tile_graph = function(tiles) {
	n   = nrow(tiles$meta_data)
	ids = tiles$meta_data$id
	adj = Matrix::sparseMatrix(
		i        = c(tiles$edges$from, tiles$edges$to),
		j        = c(tiles$edges$to,   tiles$edges$from),
		x        = 1,
		dims     = c(n, n),
		dimnames = list(ids, ids)
	)
	as(adj, "dgCMatrix")
}


#' Agglomerate initial DMT tiles into final spatial tiles
#'
#' Takes the initial tile assignment from [run_dmt()] and merges tiles using
#' two rounds of single-linkage agglomerative clustering:
#' 1. **Main merge** — tiles are merged until none exceed `max_npts` cells.
#'    Each candidate merge is scored by transcriptional similarity (GMM-based),
#'    size penalty, and change in compactness.
#' 2. **Cleanup** — remaining tiles with fewer than `min_npts` cells are
#'    absorbed into their best-scoring neighbour.
#'
#' Typical usage:
#' ```r
#' field      = compute_field(cells, mesh)
#' mesh$morse = compute_morse(field, mesh)
#' dmt        = run_dmt(mesh)
#' tiles      = make_tiles(cells, mesh, dmt, params)
#' ```
#'
#' @param cells A cells list as returned by [make_cells()], with
#'   `cells$embeddings$pca` set to an `nrow(mesh$pts)` × D matrix and
#'   (optionally) `cells$counts` set to a genes × N sparse matrix.
#' @param mesh A mesh list as returned by [make_mesh()], with `mesh$morse`
#'   populated by [compute_morse()].
#' @param dmt A dmt list as returned by [run_dmt()].
#' @param params Named list of agglomeration parameters:
#' * `alpha` — GMM sigma multiplier for transcriptional similarity scoring.
#'   `0.2` = conservative merging, `2` = liberal merging. Default `1`.
#' * `max_npts` — maximum cells per tile during the main merge. Default `50`.
#' * `min_npts` — minimum cells per tile during the cleanup pass. Default `5`.
#' * `seed` — random seed for reproducible GMM fitting. Default `1`.
#'
#' @returns A tiles list with slots:
#' * `meta_data` — data.table of tile metadata (id, X, Y, npts, shape,
#'   area, perimeter)
#' * `edges`     — data.table of inter-tile edges (from, to, coordinates,
#'   area, npts, edge_length, scores)
#' * `pcs`       — num_tiles × D matrix of per-tile mean embeddings
#' * `counts`    — genes × num_tiles sparse count matrix
#' * `adj`       — symmetric dgCMatrix tile adjacency (from [make_tile_graph()])
#' * `cell_ids`  — data.table with columns `ORIG_ID` and `tile_id` mapping
#'   each surviving cell to its final tile
#' * `params`    — list of agglomeration parameters used
#' * `meta`      — list()
#'
#' @seealso [run_dmt()], [make_tile_graph()]
#' @export
make_tiles = function(cells, mesh, dmt, params = list()) {
	alpha    = params$alpha    %||% 1
	max_npts = params$max_npts %||% 50L
	min_npts = params$min_npts %||% 5L
	seed     = params$seed     %||% 1L

	# ── Pruned counts (genes × pruned_cells) ──────────────────────────────────
	counts_all = cells$counts %||% Matrix::sparseMatrix(
		i = integer(0), j = integer(0), dims = c(0L, nrow(mesh$pts))
	)
	pruned_counts = counts_all[, mesh$pts$ORIG_ID, drop = FALSE]

	# ── Augmented dmt-like object for old-API internal functions ──────────────
	# data.table::copy ensures make_tiles never mutates the caller's dmt slots.
	#
	# The C++ functions (trace_polygons_cpp, get_agg_to_edge, …) index edges by
	# fixed column positions defined in src/shared.cpp:
	#   cols  0-14 : original mesh edge columns (from_pt … boundary)
	#   col  15    : f_prim   (IDX_F_PRIM)
	#   col  16    : f_dual   (IDX_F_DUAL)
	#   col  17    : agg_from (IDX_AGG_FROM)
	#   col  18    : agg_to   (IDX_AGG_TO)
	# run_dmt appended agg_from/agg_to directly after boundary (at 15/16), so we
	# must insert f_prim and f_dual first to push agg_from/agg_to to 17/18.
	dmt_tmp_edges = data.table::copy(mesh$edges)          # cols 0-14
	dmt_tmp_edges[, f_prim   := mesh$morse$f$edges_prim]  # col 15
	dmt_tmp_edges[, f_dual   := mesh$morse$f$edges_dual]  # col 16
	dmt_tmp_edges[, agg_from := dmt$edges$agg_from]       # col 17
	dmt_tmp_edges[, agg_to   := dmt$edges$agg_to]         # col 18

	dmt_tmp = list(
		pts       = data.table::copy(dmt$pts),
		tris      = mesh$triangles,
		edges     = dmt_tmp_edges,
		counts    = pruned_counts,
		udv_cells = list(embeddings = cells$embeddings$pca)
	)

	# ── Initialise tiles from initial DMT segmentation ────────────────────────
	aggs = dmt_init_tiles(dmt_tmp)

	# ── Round 1: merge until no tile exceeds max_npts ─────────────────────────
	aggs    = init_scores(aggs, agg_mode = 2, seed = seed, alpha = alpha, max_npts = max_npts)
	aggs    = merge_aggs(aggs, agg_mode = 2, max_npts = max_npts)
	dmt_tmp = update_dmt_aggid(dmt_tmp, aggs)
	aggs    = update_agg_shapes(dmt_tmp, aggs)

	# ── Round 2: absorb stray small tiles ─────────────────────────────────────
	aggs    = init_scores(aggs, agg_mode = 3, seed = seed, alpha = alpha, min_npts = min_npts)
	aggs    = merge_aggs(aggs, agg_mode = 3, min_npts = min_npts)
	dmt_tmp = update_dmt_aggid(dmt_tmp, aggs)
	aggs    = update_agg_shapes(dmt_tmp, aggs)
	aggs$meta_data$shape = sf::st_cast(aggs$meta_data$shape, "MULTIPOLYGON")

	# ── Build tile adjacency dgCMatrix ────────────────────────────────────────
	aggs$adj = make_tile_graph(aggs)

	# ── Per-cell tile mapping ─────────────────────────────────────────────────
	aggs$cell_ids = dmt_tmp$pts[, .(ORIG_ID, tile_id = agg_id)]

	# ── Clean up internal agglomeration state ─────────────────────────────────
	aggs$pcs_merged = NULL
	aggs$aggmap     = NULL
	aggs$d_mu       = NULL
	aggs$d_sig      = NULL

	aggs$params = list(alpha = alpha, max_npts = max_npts, min_npts = min_npts)
	aggs$meta   = list()

	aggs
}
