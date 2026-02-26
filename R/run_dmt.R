#' Assign cells to initial DMT tiles using separatrix edges
#'
#' Removes the separatrix edges stored in `mesh$morse$e_sep` from the mesh
#' graph and labels each connected component as a distinct tile.  The result
#' is used as input to [make_tiles()], which performs the agglomerative
#' merging step.
#'
#' Typical usage:
#' ```r
#' field      = compute_field(cells, mesh)
#' mesh$morse = compute_morse(field, mesh)
#' dmt        = run_dmt(mesh)
#' tiles      = make_tiles(cells, mesh, dmt, params)
#' ```
#'
#' @param mesh A mesh list as returned by [make_mesh()], with `mesh$morse`
#'   populated by [compute_morse()].
#' @param params Named list of DMT parameters (currently unused; reserved for
#'   future extensions). Default `list()`.
#'
#' @returns A dmt list with slots:
#' * `pts`   — data.table copy of `mesh$pts` with an added `agg_id` integer
#'   column giving each cell's initial tile assignment.
#' * `edges` — data.table copy of `mesh$edges` with added `agg_from` and
#'   `agg_to` integer columns giving the tile IDs of each edge's endpoints.
#' * `meta`  — list()
#'
#' @seealso [compute_morse()], [make_tiles()]
#' @export
run_dmt = function(mesh, params = list()) {
	if (is.null(mesh$morse))
		stop(paste0(
			"mesh$morse is NULL.\\n",
			"Run compute_morse() first:\\n",
			"  field      = compute_field(cells, mesh)\\n",
			"  mesh$morse = compute_morse(field, mesh)"
		))

	n_pts   = nrow(mesh$pts)
	n_edges = nrow(mesh$edges)
	e_sep   = mesh$morse$e_sep

	# ── Build reduced graph (non-separatrix edges only) ───────────────────────
	e_keep   = setdiff(seq_len(n_edges), e_sep)
	from_e   = mesh$edges$from_pt[e_keep]
	to_e     = mesh$edges$to_pt[e_keep]
	edge_vec = c(rbind(from_e, to_e))  # interleaved vertex pairs for add_edges

	g = igraph::make_empty_graph(n = n_pts, directed = FALSE)
	if (length(edge_vec) > 0) g = igraph::add_edges(g, edge_vec)

	# ── Connected components → tile IDs ──────────────────────────────────────
	agg_id = as.integer(igraph::components(g)$membership)

	# ── Attach assignments; data.table::copy prevents mutating mesh slots ─────
	pts   = data.table::copy(mesh$pts)
	edges = data.table::copy(mesh$edges)
	pts[,   agg_id   := agg_id]
	edges[, agg_from := pts$agg_id[from_pt]]
	edges[, agg_to   := pts$agg_id[to_pt]]

	list(pts = pts, edges = edges, meta = list())
}
