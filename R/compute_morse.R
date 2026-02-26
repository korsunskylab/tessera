#' Compute Discrete Morse Theory segmentation from gradient field and mesh
#'
#' Derives a scalar field from the gradient magnitude at every mesh element,
#' constructs a primal minimum spanning forest on points and a dual maximum
#' spanning forest on triangles, then traces separatrix edges that partition
#' the mesh into topologically coherent tiles.
#'
#' The caller assigns the result back to the mesh:
#' ```r
#' mesh$morse = compute_morse(field, mesh, params)
#' ```
#' The `mesh$morse` slot is then read by [run_dmt()].
#'
#' @param field A field list as returned by [compute_field()].
#' @param mesh A mesh list as returned by [make_mesh()].
#' @param params Named list of DMT parameters (currently unused; reserved for
#'   future extensions). Default `list()`.
#'
#' @returns A morse list with slots:
#' * `e_sep`  — integer vector of separatrix edge indices (1-indexed)
#' * `prim`   — primal minimum spanning forest (from [do_primary_forest()])
#' * `dual`   — dual maximum spanning forest (from [do_dual_forest()])
#' * `f`      — list of scalar field values:
#'   * `pts`        — numeric vector length N (sum of singular values at points)
#'   * `tris`       — numeric vector length F (sum of singular values at triangles)
#'   * `edges_prim` — numeric vector length E (sum of singular values at primal edges)
#'   * `edges_dual` — numeric vector length E (sum of singular values at dual edges)
#' * `params` — list() (reserved)
#' * `meta`   — list()
#'
#' @seealso [compute_field()], [run_dmt()]
#' @export
compute_morse = function(field, mesh, params = list()) {
	# ── Scalar field: sum of singular values (Frobenius norm proxy) ───────────
	f_pts        = rowSums(field$gradient$pts[,        c('len_grad', 'len_ortho')])
	f_tris       = rowSums(field$gradient$tris[,       c('len_grad', 'len_ortho')])
	f_edges_prim = rowSums(field$gradient$edges_pts[,  c('len_grad', 'len_ortho')])
	f_edges_dual = rowSums(field$gradient$edges_tris[, c('len_grad', 'len_ortho')])

	# ── Minimal dmt-like object for existing internal functions ───────────────
	# data.table::copy prevents mutation of mesh$pts / mesh$triangles / mesh$edges
	dmt_tmp = list(
		pts   = data.table::copy(mesh$pts),
		tris  = data.table::copy(mesh$triangles),
		edges = data.table::copy(mesh$edges)
	)
	dmt_tmp$pts$f        = f_pts
	dmt_tmp$tris$f       = f_tris
	dmt_tmp$edges$f_prim = f_edges_prim
	dmt_tmp$edges$f_dual = f_edges_dual

	# ── Primal minimum spanning forest on points ──────────────────────────────
	prim = do_primary_forest(dmt_tmp)

	# ── Dual maximum spanning forest on triangles ─────────────────────────────
	dual = do_dual_forest(dmt_tmp)

	# ── Separatrices ──────────────────────────────────────────────────────────
	dmt_tmp$prim = prim
	dmt_tmp$dual = dual
	e_sep = dmt_get_separatrices(dmt_tmp)

	list(
		e_sep = e_sep,
		prim  = prim,
		dual  = dual,
		f     = list(
			pts        = f_pts,
			tris       = f_tris,
			edges_prim = f_edges_prim,
			edges_dual = f_edges_dual
		),
		params = list(),
		meta   = list()
	)
}
