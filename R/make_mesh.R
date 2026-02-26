#' Construct Delaunay triangulation mesh from cells
#'
#' Builds a pruned Delaunay triangulation over cell coordinates, adds
#' degenerate exterior triangles along the boundary, and pre-computes
#' the sparse point adjacency matrix used by subsequent pipeline stages.
#'
#' Pruning has two steps: (1) edges longer than a length threshold are
#' removed along with any triangles that contain them, and (2) connected
#' components with fewer than `prune_mincells` cells are removed.
#'
#' @param cells A cells list as returned by [make_cells()].
#' @param params Named list of mesh-construction parameters:
#' * `prune_thresh_quantile` — quantile of edge length above which edges
#'   are pruned (default 0.95).
#' * `prune_mincells` — minimum cells required to keep a connected component
#'   (default 10).
#' * `prune_thresh` — absolute edge-length threshold; overrides
#'   `prune_thresh_quantile` when not `NA` (default NA).
#'
#' @returns A mesh list with slots:
#' * `pts`       — data.table of surviving cells (columns X, Y, ORIG_ID, ...)
#' * `edges`     — data.table of mesh edges (from_pt, to_pt, from_tri, to_tri, ...)
#' * `triangles` — data.table of triangles (X, Y, area, height, external)
#' * `tri_to_pt` — dgCMatrix (num_tris × num_pts): 1 where triangle uses point
#' * `adj`       — dgCMatrix (num_pts × num_pts): symmetric point adjacency
#' * `weights`   — numeric vector length num_edges (all 1, reserved for future use)
#' * `morse`     — NULL until populated by [compute_morse()]
#' * `meta`      — list()
#'
#' @seealso [make_cells()], [compute_field()]
#' @export
make_mesh = function(cells, params = list()) {
	prune_thresh_quantile = params$prune_thresh_quantile %||% 0.95
	prune_mincells        = params$prune_mincells        %||% 10
	prune_thresh          = params$prune_thresh          %||% NA

	# Placeholder counts when cells carries none (prune_graph needs a counts slot)
	counts = cells$counts %||% Matrix::sparseMatrix(
		i = integer(0), j = integer(0),
		dims = c(0L, nrow(cells$coords))
	)

	# Build triangulation, prune long edges / small components, add boundary tris.
	# init_data / prune_graph / add_exterior_triangles are the canonical
	# implementations; they will be inlined here when utils_initdata.R is
	# retired in Block 7.
	data = init_data(cells$coords[, 1], cells$coords[, 2], counts)
	data = prune_graph(data, prune_thresh_quantile, prune_mincells, prune_thresh)
	data = add_exterior_triangles(data)

	# Build sparse point adjacency matrix (num_pts x num_pts, dgCMatrix).
	# Diagonal is left as 0; compute_field sets it to 1 when smoothing.
	n_pts    = nrow(data$pts)
	edge_vec = c(t(as.matrix(data$edges[, .(from_pt, to_pt)])))
	g        = igraph::make_empty_graph(n = n_pts, directed = FALSE)
	if (length(edge_vec) > 0) g = igraph::add_edges(g, edge_vec)
	adj = as(igraph::as_adjacency_matrix(g, sparse = TRUE), "dgCMatrix")

	list(
		pts       = data$pts,
		edges     = data$edges,
		triangles = data$tris,
		tri_to_pt = data$tri_to_pt,
		adj       = adj,
		weights   = rep(1, nrow(data$edges)),
		morse     = NULL,
		meta      = list()
	)
}
