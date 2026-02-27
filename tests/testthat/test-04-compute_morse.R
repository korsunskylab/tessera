# ── Structure tests ────────────────────────────────────────────────────────────

test_that("compute_morse returns list with correct top-level slots", {
	field = compute_field(.test_cells_field, .test_mesh)
	morse = compute_morse(field, .test_mesh)
	expect_type(morse, "list")
	expect_named(morse, c("e_sep", "prim", "dual", "f", "params", "meta"))
})

test_that("compute_morse f has correct sub-slots", {
	field = compute_field(.test_cells_field, .test_mesh)
	morse = compute_morse(field, .test_mesh)
	expect_named(morse$f, c("pts", "tris", "edges_prim", "edges_dual"))
})

test_that("compute_morse prim has correct sub-slots", {
	field = compute_field(.test_cells_field, .test_mesh)
	morse = compute_morse(field, .test_mesh)
	expect_named(morse$prim, c("edges", "saddles", "labels", "minima", "parent", "parent_edge"))
})

test_that("compute_morse dual has correct sub-slots", {
	field = compute_field(.test_cells_field, .test_mesh)
	morse = compute_morse(field, .test_mesh)
	expect_named(morse$dual, c("edges", "saddles", "labels", "maxima", "parent", "parent_edge"))
})

test_that("compute_morse f vectors have correct lengths", {
	field  = compute_field(.test_cells_field, .test_mesh)
	morse  = compute_morse(field, .test_mesh)
	n_pts  = nrow(.test_mesh$pts)
	n_tris = nrow(.test_mesh$triangles)
	n_edges = nrow(.test_mesh$edges)

	expect_length(morse$f$pts,        n_pts)
	expect_length(morse$f$tris,       n_tris)
	expect_length(morse$f$edges_prim, n_edges)
	expect_length(morse$f$edges_dual, n_edges)
})

test_that("compute_morse f values are non-negative", {
	field = compute_field(.test_cells_field, .test_mesh)
	morse = compute_morse(field, .test_mesh)
	expect_true(all(morse$f$pts        >= 0))
	expect_true(all(morse$f$tris       >= 0))
	expect_true(all(morse$f$edges_prim >= 0))
	expect_true(all(morse$f$edges_dual >= 0))
})

test_that("compute_morse e_sep is an integer vector of valid edge indices", {
	field   = compute_field(.test_cells_field, .test_mesh)
	morse   = compute_morse(field, .test_mesh)
	n_edges = nrow(.test_mesh$edges)
	expect_type(morse$e_sep, "integer")
	expect_true(length(morse$e_sep) > 0)
	expect_true(all(morse$e_sep >= 1L & morse$e_sep <= n_edges))
})

test_that("compute_morse prim and dual labels have correct lengths", {
	field  = compute_field(.test_cells_field, .test_mesh)
	morse  = compute_morse(field, .test_mesh)
	n_pts  = nrow(.test_mesh$pts)
	n_tris = nrow(.test_mesh$triangles)
	expect_length(morse$prim$labels, n_pts)
	expect_length(morse$dual$labels, n_tris)
})

test_that("compute_morse params and meta are empty lists", {
	field = compute_field(.test_cells_field, .test_mesh)
	morse = compute_morse(field, .test_mesh)
	expect_type(morse$params, "list")
	expect_length(morse$params, 0L)
	expect_type(morse$meta, "list")
	expect_length(morse$meta, 0L)
})

test_that("compute_morse does not mutate mesh$pts", {
	field      = compute_field(.test_cells_field, .test_mesh)
	pts_before = names(.test_mesh$pts)
	compute_morse(field, .test_mesh)
	expect_equal(names(.test_mesh$pts), pts_before)
	expect_false("f" %in% names(.test_mesh$pts))
})

