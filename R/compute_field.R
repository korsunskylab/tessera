#' Compute spatial gradient field from cell embeddings and mesh geometry
#'
#' Estimates a spatial gradient at every mesh element (points, triangles,
#' edges) from cell embeddings, applies optional bilateral smoothing, and
#' returns SVD-compressed gradient representations used by [compute_morse()].
#'
#' **Important:** `cells$embeddings$pca` must contain exactly `nrow(mesh$pts)`
#' rows — one per cell that survived pruning in [make_mesh()]. Compute
#' embeddings on the pruned cells after [make_mesh()]:
#' ```r
#' mesh             = make_mesh(cells, params)
#' pca              = do_pca(cells$counts[, mesh$pts$ORIG_ID], npcs = 10)
#' cells$embeddings = list(pca = pca$embeddings)
#' field            = compute_field(cells, mesh, params)
#' ```
#'
#' @param cells A cells list as returned by [make_cells()], with
#'   `cells$embeddings$pca` set to an `nrow(mesh$pts)` × D matrix.
#' @param mesh A mesh list as returned by [make_mesh()].
#' @param params Named list of gradient parameters:
#' * `smooth_distance` — distance weighting for bilateral smoothing.
#'   One of `'none'`, `'euclidean'`, `'projected'`, `'constant'`.
#'   Default `'projected'`.
#' * `smooth_similarity` — similarity weighting for bilateral smoothing.
#'   One of `'none'`, `'euclidean'`, `'projected'`, `'constant'`.
#'   Default `'projected'`.
#' * `smooth_iter` — number of smoothing iterations. Default `1`.
#'
#' @returns A field list with slots:
#' * `gradient` — list with:
#'   * `pts`        — N × 6 matrix of SVD-compressed point gradients
#'   * `tris`       — F × 6 matrix of SVD-compressed triangle gradients
#'   * `edges_pts`  — E × 6 matrix of SVD-compressed primal edge gradients
#'   * `edges_tris` — E × 6 matrix of SVD-compressed dual edge gradients
#'   * `params`     — list of smoothing parameters used
#'   * `meta`       — list()
#'
#'   All SVD matrices have columns:
#'   `dx_grad`, `dy_grad`, `dx_ortho`, `dy_ortho`, `len_grad`, `len_ortho`.
#' * `meta` — list()
#'
#' @seealso [make_mesh()], [compute_morse()]
#' @export
compute_field = function(cells, mesh, params = list()) {
	smooth_distance  = params$smooth_distance   %||% 'projected'
	smooth_similarity = params$smooth_similarity %||% 'projected'
	smooth_iter      = params$smooth_iter        %||% 1L

	embeddings = cells$embeddings$pca
	if (is.null(embeddings))
		stop(paste0(
			"cells$embeddings$pca is NULL.\n",
			"Compute PCA on pruned cells after make_mesh():\n",
			"  pca              = do_pca(cells$counts[, mesh$pts$ORIG_ID], npcs)\n",
			"  cells$embeddings = list(pca = pca$embeddings)"
		))
	if (nrow(embeddings) != nrow(mesh$pts))
		stop(paste0(
			"cells$embeddings$pca has ", nrow(embeddings),
			" rows but mesh has ", nrow(mesh$pts), " points.\n",
			"Recompute PCA on pruned cells:\n",
			"  pca              = do_pca(cells$counts[, mesh$pts$ORIG_ID], npcs)\n",
			"  cells$embeddings = list(pca = pca$embeddings)"
		))

	coords = as.matrix(mesh$pts[, .(X, Y)])

	# ── Point gradients ───────────────────────────────────────────────────────
	# estimate_field sets diag(adj) = 0 internally on a copy; mesh$adj unmodified
	grad_pts = estimate_field(coords, mesh$adj, embeddings)

	# ── Optional bilateral smoothing ──────────────────────────────────────────
	# smooth_field sets diag(adj) = 1 internally on a copy; mesh$adj unmodified
	if (smooth_distance != 'none' && smooth_similarity != 'none') {
		for (i in seq_len(smooth_iter)) {
			grad_pts = smooth_field(
				coords,
				field      = grad_pts,
				adj        = mesh$adj,
				include_self = TRUE,
				distance   = smooth_distance,
				similarity = smooth_similarity
			)
		}
	}

	# ── Triangle gradients (weighted average of vertex gradients) ─────────────
	n_tris = nrow(mesh$triangles)
	D      = dim(grad_pts)[2]
	grad_tris        = array(dim = c(2L, D, n_tris))
	grad_tris[1L, ,] = as.matrix(grad_pts[1L, ,] %*% Matrix::t(mesh$tri_to_pt))
	grad_tris[2L, ,] = as.matrix(grad_pts[2L, ,] %*% Matrix::t(mesh$tri_to_pt))
	grad_tris        = sweep(grad_tris, 3L, Matrix::rowSums(mesh$tri_to_pt), '/')

	# ── Primal edge gradients (average of two endpoint gradients) ─────────────
	n_pts = nrow(mesh$pts)
	y1    = factor(mesh$edges$from_pt, levels = seq_len(n_pts))
	y2    = factor(mesh$edges$to_pt,   levels = seq_len(n_pts))
	adj_e_to_pts     = Matrix::sparse.model.matrix(~0 + y1) +
	                   Matrix::sparse.model.matrix(~0 + y2)
	grad_edges_pts        = array(dim = c(2L, D, nrow(mesh$edges)))
	grad_edges_pts[1L, ,] = as.matrix(grad_pts[1L, ,] %*% Matrix::t(adj_e_to_pts))
	grad_edges_pts[2L, ,] = as.matrix(grad_pts[2L, ,] %*% Matrix::t(adj_e_to_pts))
	grad_edges_pts        = grad_edges_pts / 2

	# ── Dual edge gradients (average of two adjacent triangle gradients) ───────
	y1    = factor(mesh$edges$from_tri, levels = seq_len(n_tris))
	y2    = factor(mesh$edges$to_tri,   levels = seq_len(n_tris))
	adj_e_to_tris     = Matrix::sparse.model.matrix(~0 + y1) +
	                    Matrix::sparse.model.matrix(~0 + y2)
	grad_edges_tris        = array(dim = c(2L, D, nrow(mesh$edges)))
	grad_edges_tris[1L, ,] = as.matrix(grad_tris[1L, ,] %*% Matrix::t(adj_e_to_tris))
	grad_edges_tris[2L, ,] = as.matrix(grad_tris[2L, ,] %*% Matrix::t(adj_e_to_tris))
	grad_edges_tris        = grad_edges_tris / 2

	# ── SVD compression at each mesh location ─────────────────────────────────
	# compress_field_cpp returns N × 6: dx_grad dy_grad dx_ortho dy_ortho len_grad len_ortho
	svd_cols = c('dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho', 'len_grad', 'len_ortho')

	pts_svd        = compress_field_cpp(grad_pts)       ; colnames(pts_svd)        = svd_cols
	tris_svd       = compress_field_cpp(grad_tris)      ; colnames(tris_svd)       = svd_cols
	edges_pts_svd  = compress_field_cpp(grad_edges_pts) ; colnames(edges_pts_svd)  = svd_cols
	edges_tris_svd = compress_field_cpp(grad_edges_tris); colnames(edges_tris_svd) = svd_cols

	list(
		gradient = list(
			pts        = pts_svd,
			tris       = tris_svd,
			edges_pts  = edges_pts_svd,
			edges_tris = edges_tris_svd,
			params     = list(
				smooth_distance   = smooth_distance,
				smooth_similarity = smooth_similarity,
				smooth_iter       = smooth_iter
			),
			meta = list()
		),
		meta = list()
	)
}
