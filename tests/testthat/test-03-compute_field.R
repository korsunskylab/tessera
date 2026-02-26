# ── Structure tests ────────────────────────────────────────────────────────────

test_that("compute_field returns list with correct top-level slots", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_type(field, "list")
	expect_named(field, c("gradient", "meta"))
})

test_that("compute_field gradient has correct sub-slots", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_named(field$gradient,
		c("pts", "tris", "edges_pts", "edges_tris", "params", "meta"))
})

test_that("compute_field gradient matrices have correct dimensions", {
	field  = compute_field(.test_cells_field, .test_mesh)
	n_pts  = nrow(.test_mesh$pts)
	n_tris = nrow(.test_mesh$triangles)
	n_edges = nrow(.test_mesh$edges)

	expect_equal(dim(field$gradient$pts),        c(n_pts,   6L))
	expect_equal(dim(field$gradient$tris),       c(n_tris,  6L))
	expect_equal(dim(field$gradient$edges_pts),  c(n_edges, 6L))
	expect_equal(dim(field$gradient$edges_tris), c(n_edges, 6L))
})

test_that("compute_field gradient matrices have correct column names", {
	field    = compute_field(.test_cells_field, .test_mesh)
	expected = c('dx_grad', 'dy_grad', 'dx_ortho', 'dy_ortho', 'len_grad', 'len_ortho')
	expect_equal(colnames(field$gradient$pts),        expected)
	expect_equal(colnames(field$gradient$tris),       expected)
	expect_equal(colnames(field$gradient$edges_pts),  expected)
	expect_equal(colnames(field$gradient$edges_tris), expected)
})

test_that("compute_field gradient singular values are non-negative", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_true(all(field$gradient$pts[, 'len_grad']  >= 0))
	expect_true(all(field$gradient$pts[, 'len_ortho'] >= 0))
	expect_true(all(field$gradient$tris[, 'len_grad'] >= 0))
})

test_that("compute_field gradient unit vectors are unit length", {
	field = compute_field(.test_cells_field, .test_mesh)
	# dx_grad^2 + dy_grad^2 should equal 1 (up to tolerance) for non-zero gradients
	len_sq = field$gradient$pts[, 'dx_grad']^2 + field$gradient$pts[, 'dy_grad']^2
	nonzero = field$gradient$pts[, 'len_grad'] > 1e-10
	expect_true(all(abs(len_sq[nonzero] - 1) < 1e-6))
})

test_that("compute_field params slot records smoothing settings", {
	params = list(smooth_distance='projected', smooth_similarity='projected', smooth_iter=2L)
	field  = compute_field(.test_cells_field, .test_mesh, params)
	expect_equal(field$gradient$params$smooth_distance,   'projected')
	expect_equal(field$gradient$params$smooth_similarity, 'projected')
	expect_equal(field$gradient$params$smooth_iter,       2L)
})

test_that("compute_field meta slots are empty lists", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_type(field$meta, "list")
	expect_length(field$meta, 0L)
	expect_type(field$gradient$meta, "list")
	expect_length(field$gradient$meta, 0L)
})

test_that("compute_field validates embeddings row count", {
	bad_cells            = .test_cells_field
	bad_cells$embeddings = list(pca = .test_pca_pruned$embeddings[1:10, ])
	expect_error(compute_field(bad_cells, .test_mesh), "rows but mesh has")
})

test_that("compute_field validates embeddings are not NULL", {
	bad_cells            = .test_cells_field
	bad_cells$embeddings = list(pca = NULL)
	expect_error(compute_field(bad_cells, .test_mesh), "NULL")
})

test_that("compute_field works with no smoothing", {
	field = compute_field(.test_cells_field, .test_mesh,
		params = list(smooth_distance = 'none', smooth_similarity = 'none'))
	expect_equal(dim(field$gradient$pts), c(nrow(.test_mesh$pts), 6L))
})

# ── Numerical equivalence tests ───────────────────────────────────────────────

test_that("compute_field pts_svd matches old compress_gradients_svd output", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_equal(field$gradient$pts, .test_field_old$pts_svd, tolerance = 1e-10)
})

test_that("compute_field tris_svd matches old compress_gradients_svd output", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_equal(field$gradient$tris, .test_field_old$tris_svd, tolerance = 1e-10)
})

test_that("compute_field edges_pts_svd matches old compress_gradients_svd output", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_equal(field$gradient$edges_pts, .test_field_old$edges_pts_svd, tolerance = 1e-10)
})

test_that("compute_field edges_tris_svd matches old compress_gradients_svd output", {
	field = compute_field(.test_cells_field, .test_mesh)
	expect_equal(field$gradient$edges_tris, .test_field_old$edges_tris_svd, tolerance = 1e-10)
})
