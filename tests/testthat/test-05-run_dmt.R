# ── Structure tests ────────────────────────────────────────────────────────────

test_that("run_dmt returns list with correct top-level slots", {
	dmt = run_dmt(.test_mesh_morse)
	expect_type(dmt, "list")
	expect_named(dmt, c("pts", "edges", "meta"))
})

test_that("run_dmt pts has agg_id column", {
	dmt = run_dmt(.test_mesh_morse)
	expect_true("agg_id" %in% names(dmt$pts))
})

test_that("run_dmt edges has agg_from and agg_to columns", {
	dmt = run_dmt(.test_mesh_morse)
	expect_true("agg_from" %in% names(dmt$edges))
	expect_true("agg_to"   %in% names(dmt$edges))
})

test_that("run_dmt meta is an empty list", {
	dmt = run_dmt(.test_mesh_morse)
	expect_type(dmt$meta, "list")
	expect_length(dmt$meta, 0L)
})

test_that("run_dmt pts has same row count as mesh$pts", {
	dmt = run_dmt(.test_mesh_morse)
	expect_equal(nrow(dmt$pts), nrow(.test_mesh_morse$pts))
})

test_that("run_dmt edges has same row count as mesh$edges", {
	dmt = run_dmt(.test_mesh_morse)
	expect_equal(nrow(dmt$edges), nrow(.test_mesh_morse$edges))
})

# ── Property tests ─────────────────────────────────────────────────────────────

test_that("run_dmt agg_id is a positive integer vector", {
	dmt = run_dmt(.test_mesh_morse)
	expect_type(dmt$pts$agg_id, "integer")
	expect_true(all(dmt$pts$agg_id >= 1L))
})

test_that("run_dmt produces at least 2 tiles", {
	dmt = run_dmt(.test_mesh_morse)
	expect_gte(length(unique(dmt$pts$agg_id)), 2L)
})

test_that("run_dmt agg_from and agg_to are consistent with pts agg_id", {
	dmt = run_dmt(.test_mesh_morse)
	expect_equal(dmt$edges$agg_from, dmt$pts$agg_id[dmt$edges$from_pt])
	expect_equal(dmt$edges$agg_to,   dmt$pts$agg_id[dmt$edges$to_pt])
})

test_that("run_dmt non-separatrix edges connect cells in the same tile", {
	dmt   = run_dmt(.test_mesh_morse)
	e_sep = .test_mesh_morse$morse$e_sep
	n_edges = nrow(dmt$edges)
	e_keep  = setdiff(seq_len(n_edges), e_sep)
	expect_true(all(dmt$edges$agg_from[e_keep] == dmt$edges$agg_to[e_keep]))
})

test_that("run_dmt some separatrix edges separate tiles", {
	dmt   = run_dmt(.test_mesh_morse)
	e_sep = .test_mesh_morse$morse$e_sep
	expect_true(any(dmt$edges$agg_from[e_sep] != dmt$edges$agg_to[e_sep]))
})

test_that("run_dmt does not mutate mesh$pts", {
	pts_cols_before = names(.test_mesh_morse$pts)
	run_dmt(.test_mesh_morse)
	expect_equal(names(.test_mesh_morse$pts), pts_cols_before)
	expect_false("agg_id" %in% names(.test_mesh_morse$pts))
})

test_that("run_dmt does not mutate mesh$edges", {
	edge_cols_before = names(.test_mesh_morse$edges)
	run_dmt(.test_mesh_morse)
	expect_equal(names(.test_mesh_morse$edges), edge_cols_before)
	expect_false("agg_from" %in% names(.test_mesh_morse$edges))
})

test_that("run_dmt errors when mesh$morse is NULL", {
	bad_mesh = .test_mesh_morse
	bad_mesh$morse = NULL
	expect_error(run_dmt(bad_mesh), "mesh\\$morse is NULL")
})

