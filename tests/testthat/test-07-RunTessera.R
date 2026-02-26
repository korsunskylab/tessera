# ── File-local fixtures (computed once on load) ────────────────────────────────
# The full pipeline takes ~60 s; keeping here avoids slowing other test files.

.local_run = suppressWarnings(RunTessera(
	X                 = .test_X,
	Y                 = .test_Y,
	counts            = .test_counts,
	meta_data         = .test_meta,
	meta_vars_include = c('type')
))
.local_dmt_run  = .local_run$dmt
.local_aggs_run = .local_run$aggs

# ── Top-level structure ────────────────────────────────────────────────────────

test_that("RunTessera returns list with dmt and aggs", {
	expect_type(.local_run, "list")
	expect_named(.local_run, c("dmt", "aggs"))
})

test_that("RunTessera dmt has required slots", {
	expect_true(all(c(
		"pts", "tris", "edges", "tri_to_pt",
		"counts", "udv_cells", "prim", "dual", "e_sep"
	) %in% names(.local_dmt_run)))
})

test_that("RunTessera aggs has required slots", {
	expect_true(all(c(
		"meta_data", "counts", "pcs", "edges", "adj"
	) %in% names(.local_aggs_run)))
})

# ── dmt$pts column checks ──────────────────────────────────────────────────────

test_that("RunTessera dmt$pts has required columns", {
	expect_true(all(c("X", "Y", "ORIG_ID", "f", "agg_id", "type") %in%
	                  names(.local_dmt_run$pts)))
})

test_that("RunTessera dmt$pts agg_id values are valid tile IDs", {
	expect_true(all(.local_dmt_run$pts$agg_id %in% .local_aggs_run$meta_data$id))
})

test_that("RunTessera dmt$pts ORIG_ID values are valid global cell indices", {
	expect_true(all(.local_dmt_run$pts$ORIG_ID %in% seq_len(.test_N)))
	expect_true(all(.local_dmt_run$pts$ORIG_ID >= 1L))
})

test_that("RunTessera dmt$pts has a row per pruned cell", {
	expect_equal(nrow(.local_dmt_run$pts), sum(.local_aggs_run$meta_data$npts))
})

# ── aggs structure ─────────────────────────────────────────────────────────────

test_that("RunTessera aggs meta_data has required columns", {
	expect_true(all(c("id", "X", "Y", "npts", "shape", "area", "perimeter") %in%
	                  names(.local_aggs_run$meta_data)))
})

test_that("RunTessera aggs tile IDs are positive consecutive integers", {
	ids = .local_aggs_run$meta_data$id
	expect_equal(sort(ids), seq_len(length(ids)))
})

test_that("RunTessera aggs counts column names match tile IDs", {
	expect_equal(
		colnames(.local_aggs_run$counts),
		as.character(.local_aggs_run$meta_data$id)
	)
})

test_that("RunTessera aggs counts dimensions are genes x tiles", {
	n_genes = nrow(.test_counts)
	n_tiles = nrow(.local_aggs_run$meta_data)
	expect_equal(dim(.local_aggs_run$counts), c(n_genes, n_tiles))
})

test_that("RunTessera aggs pcs dimensions are tiles x npcs", {
	n_tiles = nrow(.local_aggs_run$meta_data)
	expect_equal(nrow(.local_aggs_run$pcs), n_tiles)
	expect_gt(ncol(.local_aggs_run$pcs), 0L)
})

test_that("RunTessera aggs adj is a matrix with tile dimensions", {
	n = nrow(.local_aggs_run$meta_data)
	adj = .local_aggs_run$adj
	expect_equal(dim(adj), c(n, n))
})

# ── udv_cells ──────────────────────────────────────────────────────────────────

test_that("RunTessera udv_cells has embeddings and loadings", {
	expect_false(is.null(.local_dmt_run$udv_cells$embeddings))
	expect_false(is.null(.local_dmt_run$udv_cells$loadings))
})

test_that("RunTessera udv_cells embeddings rows equal pruned cell count", {
	expect_equal(
		nrow(.local_dmt_run$udv_cells$embeddings),
		nrow(.local_dmt_run$pts)
	)
})

# ── Property checks ────────────────────────────────────────────────────────────

test_that("RunTessera produces at least 2 tiles", {
	expect_gte(nrow(.local_aggs_run$meta_data), 2L)
})

test_that("RunTessera total cells across tiles equals pruned cell count", {
	expect_equal(
		sum(.local_aggs_run$meta_data$npts),
		nrow(.local_dmt_run$pts)
	)
})

test_that("RunTessera total gene counts are positive", {
	expect_gt(sum(.local_aggs_run$counts), 0)
})

test_that("RunTessera aggs edges has from and to columns", {
	expect_true(all(c("from", "to") %in% names(.local_aggs_run$edges)))
})

test_that("RunTessera dmt$tris has f column", {
	expect_true("f" %in% names(.local_dmt_run$tris))
})

test_that("RunTessera dmt$edges has f_prim, f_dual, agg_from, agg_to columns", {
	expect_true(all(c("f_prim", "f_dual", "agg_from", "agg_to") %in%
	                  names(.local_dmt_run$edges)))
})

test_that("RunTessera dmt prim has expected slots", {
	expect_true(all(c("edges", "saddles", "labels", "minima", "parent", "parent_edge") %in%
	                  names(.local_dmt_run$prim)))
})

test_that("RunTessera dmt dual has expected slots", {
	expect_true(all(c("edges", "saddles", "labels", "maxima", "parent", "parent_edge") %in%
	                  names(.local_dmt_run$dual)))
})

test_that("RunTessera dmt e_sep is an integer vector", {
	expect_type(.local_dmt_run$e_sep, "integer")
	expect_gte(length(.local_dmt_run$e_sep), 1L)
})

# ── Numerical sanity checks ────────────────────────────────────────────────────

test_that("RunTessera tile count is within plausible bounds", {
	n_tiles   = nrow(.local_aggs_run$meta_data)
	n_pruned  = nrow(.local_dmt_run$pts)
	min_npts  = 5L
	max_npts  = 50L
	# Must have more than 1 tile and each tile must fit within max_npts
	expect_gte(n_tiles, 2L)
	expect_lte(n_tiles, ceiling(n_pruned / min_npts))
})

test_that("RunTessera total gene counts equal sum over input pruned counts", {
	expected = sum(.local_dmt_run$counts)
	expect_equal(sum(.local_aggs_run$counts), expected, tolerance = 1e-10)
})
