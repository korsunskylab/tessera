# ── File-local fixtures (computed once on load, not in shared helper) ─────────
# Agglomeration takes ~30-60 s; keeping these here avoids slowing every other
# test file.

.local_dmt   = run_dmt(.test_mesh_morse)
.local_tiles = make_tiles(.test_cells_field, .test_mesh_morse, .local_dmt)

# Old-pipeline reference: build dmt-like object from existing fixtures then run
# the same agglomeration calls so we can compare structural outputs.
.local_tiles_old = local({
	dmt = list(
		pts = local({
			p        = data.table::copy(.test_morse_old$pts)
			p$agg_id = as.integer(.test_dmt_old$agg_id)
			p
		}),
		tris  = .test_mesh_old$tris,
		edges = local({
			e          = data.table::copy(.test_morse_old$edges)
			e$agg_from = as.integer(.test_dmt_old$agg_id)[e$from_pt]
			e$agg_to   = as.integer(.test_dmt_old$agg_id)[e$to_pt]
			e
		}),
		counts    = .test_counts[, .test_mesh$pts$ORIG_ID, drop = FALSE],
		udv_cells = list(embeddings = .test_pca_pruned$embeddings)
	)
	aggs    = dmt_init_tiles(dmt)
	aggs    = init_scores(aggs, agg_mode = 2, alpha = 1, max_npts = 50L)
	aggs    = merge_aggs(aggs,  agg_mode = 2, max_npts = 50L)
	dmt     = update_dmt_aggid(dmt, aggs)
	aggs    = update_agg_shapes(dmt, aggs)
	aggs    = init_scores(aggs, agg_mode = 3, alpha = 1, min_npts = 5L)
	aggs    = merge_aggs(aggs,  agg_mode = 3, min_npts = 5L)
	dmt     = update_dmt_aggid(dmt, aggs)
	aggs    = update_agg_shapes(dmt, aggs)
	aggs$meta_data$shape = sf::st_cast(aggs$meta_data$shape, "MULTIPOLYGON")
	aggs
})

# ── make_tile_graph tests ──────────────────────────────────────────────────────

test_that("make_tile_graph returns a dgCMatrix", {
	adj = make_tile_graph(.local_tiles)
	expect_s4_class(adj, "dgCMatrix")
})

test_that("make_tile_graph dimensions match number of tiles", {
	adj     = make_tile_graph(.local_tiles)
	n_tiles = nrow(.local_tiles$meta_data)
	expect_equal(dim(adj), c(n_tiles, n_tiles))
})

test_that("make_tile_graph is symmetric", {
	adj = make_tile_graph(.local_tiles)
	expect_equal(adj, Matrix::t(adj))
})

test_that("make_tile_graph diagonal is zero", {
	adj = make_tile_graph(.local_tiles)
	expect_true(all(Matrix::diag(adj) == 0))
})

# ── make_tiles structure tests ─────────────────────────────────────────────────

test_that("make_tiles returns list with correct top-level slots", {
	expect_type(.local_tiles, "list")
	expect_true(all(c("meta_data", "edges", "pcs", "counts", "adj",
	                  "cell_ids", "params", "meta") %in% names(.local_tiles)))
})

test_that("make_tiles meta_data has required columns", {
	expect_true(all(c("id", "X", "Y", "npts", "shape", "area", "perimeter") %in%
	                  names(.local_tiles$meta_data)))
})

test_that("make_tiles adj is a dgCMatrix", {
	expect_s4_class(.local_tiles$adj, "dgCMatrix")
})

test_that("make_tiles adj dimensions match number of tiles", {
	n = nrow(.local_tiles$meta_data)
	expect_equal(dim(.local_tiles$adj), c(n, n))
})

test_that("make_tiles cell_ids has ORIG_ID and tile_id columns", {
	expect_named(.local_tiles$cell_ids, c("ORIG_ID", "tile_id"))
})

test_that("make_tiles cell_ids has one row per pruned cell", {
	expect_equal(nrow(.local_tiles$cell_ids), nrow(.test_mesh_morse$pts))
})

test_that("make_tiles params records alpha, max_npts, min_npts", {
	expect_equal(.local_tiles$params$alpha,    1)
	expect_equal(.local_tiles$params$max_npts, 50L)
	expect_equal(.local_tiles$params$min_npts, 5L)
})

test_that("make_tiles meta is an empty list", {
	expect_type(.local_tiles$meta, "list")
	expect_length(.local_tiles$meta, 0L)
})

# ── make_tiles property tests ──────────────────────────────────────────────────

test_that("make_tiles produces at least 2 tiles", {
	expect_gte(nrow(.local_tiles$meta_data), 2L)
})

test_that("make_tiles total cells equals number of pruned mesh points", {
	expect_equal(sum(.local_tiles$meta_data$npts), nrow(.test_mesh_morse$pts))
})

test_that("make_tiles all tiles have at least min_npts cells", {
	expect_true(all(.local_tiles$meta_data$npts >= .local_tiles$params$min_npts))
})

test_that("make_tiles no tile exceeds max_npts cells", {
	expect_true(all(.local_tiles$meta_data$npts <= .local_tiles$params$max_npts))
})

test_that("make_tiles pcs dimensions are tiles x npcs", {
	n_tiles = nrow(.local_tiles$meta_data)
	npcs    = ncol(.test_pca_pruned$embeddings)
	expect_equal(dim(.local_tiles$pcs), c(n_tiles, npcs))
})

test_that("make_tiles counts dimensions are genes x tiles", {
	n_genes = nrow(.test_counts)
	n_tiles = nrow(.local_tiles$meta_data)
	expect_equal(dim(.local_tiles$counts), c(n_genes, n_tiles))
})

test_that("make_tiles tile IDs are positive consecutive integers", {
	ids = .local_tiles$meta_data$id
	expect_equal(sort(ids), seq_len(length(ids)))
})

test_that("make_tiles cell_ids tile_id values are valid tile IDs", {
	valid_ids = .local_tiles$meta_data$id
	expect_true(all(.local_tiles$cell_ids$tile_id %in% valid_ids))
})

test_that("make_tiles does not mutate dmt$pts", {
	pts_cols_before = names(.local_dmt$pts)
	agg_id_before   = .local_dmt$pts$agg_id
	# re-run make_tiles with the same dmt
	make_tiles(.test_cells_field, .test_mesh_morse, .local_dmt)
	expect_equal(names(.local_dmt$pts),      pts_cols_before)
	expect_equal(.local_dmt$pts$agg_id,      agg_id_before)
})

# ── Numerical equivalence tests ───────────────────────────────────────────────

test_that("make_tiles produces same number of tiles as old pipeline", {
	expect_equal(nrow(.local_tiles$meta_data), nrow(.local_tiles_old$meta_data))
})

test_that("make_tiles total gene counts match old pipeline", {
	expect_equal(
		sum(.local_tiles$counts),
		sum(.local_tiles_old$counts),
		tolerance = 1e-10
	)
})
