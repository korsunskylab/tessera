test_that("make_cells returns a list with correct slot names", {
	cells = make_cells(.test_coords, .test_embeddings, .test_counts, .test_meta)
	expect_type(cells, "list")
	expect_named(cells, c("coords", "embeddings", "counts", "meta", "params"))
})

test_that("make_cells coords is N x 2 matrix with named columns", {
	cells = make_cells(.test_coords, .test_embeddings, .test_counts, .test_meta)
	expect_true(is.matrix(cells$coords))
	expect_equal(nrow(cells$coords), .test_N)
	expect_equal(ncol(cells$coords), 2)
	expect_equal(colnames(cells$coords), c("x", "y"))
	expect_equal(unname(cells$coords[, "x"]), .test_X)
	expect_equal(unname(cells$coords[, "y"]), .test_Y)
})

test_that("make_cells counts is dgCMatrix with N columns", {
	cells = make_cells(.test_coords, .test_embeddings, .test_counts, .test_meta)
	expect_s4_class(cells$counts, "dgCMatrix")
	expect_equal(ncol(cells$counts), .test_N)
	expect_equal(nrow(cells$counts), nrow(.test_counts))
})

test_that("make_cells meta is data.table with N rows", {
	cells = make_cells(.test_coords, .test_embeddings, .test_counts, .test_meta)
	expect_s3_class(cells$meta, "data.table")
	expect_equal(nrow(cells$meta), .test_N)
})

test_that("make_cells embeddings match do_pca output numerically", {
	cells = make_cells(.test_coords, .test_embeddings, .test_counts, .test_meta)
	expect_true(is.list(cells$embeddings))
	expect_named(cells$embeddings, "pca")
	expect_equal(cells$embeddings$pca, .test_pca$embeddings)
})

test_that("make_cells works without optional counts, meta, params", {
	cells = make_cells(.test_coords, .test_embeddings)
	expect_null(cells$counts)
	expect_null(cells$meta)
	expect_type(cells$params, "list")
	expect_length(cells$params, 0)
})

test_that("make_cells coerces data.frame meta to data.table", {
	cells = make_cells(.test_coords, .test_embeddings, meta = as.data.frame(.test_meta))
	expect_s3_class(cells$meta, "data.table")
})

test_that("make_cells validates coords shape", {
	expect_error(make_cells(matrix(1:6, ncol = 3), .test_embeddings), "N x 2")
	expect_error(make_cells(as.data.frame(.test_coords), .test_embeddings), "N x 2")
})

test_that("make_cells validates embeddings row count", {
	bad = list(pca = matrix(1:20, nrow = 5))  # 5 rows, not N
	expect_error(make_cells(.test_coords, bad), "N rows")
})

test_that("make_cells validates embeddings is named list", {
	expect_error(make_cells(.test_coords, list(.test_pca$embeddings)), "named list")
})

test_that("make_cells validates counts column count", {
	bad_counts = .test_counts[, 1:10]  # 10 cells, not N
	expect_error(make_cells(.test_coords, .test_embeddings, counts = bad_counts), "N columns")
})

test_that("make_cells validates meta row count", {
	bad_meta = .test_meta[1:10, ]
	expect_error(make_cells(.test_coords, .test_embeddings, meta = bad_meta), "N rows")
})
