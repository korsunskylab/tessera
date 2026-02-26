# ── Structure tests ────────────────────────────────────────────────────────────

test_that("make_mesh returns list with correct slot names", {
	mesh = make_mesh(.test_cells)
	expect_type(mesh, "list")
	expect_named(mesh, c("pts", "edges", "triangles", "tri_to_pt", "adj", "weights", "morse", "meta"))
})

test_that("make_mesh pts is data.table with X, Y, ORIG_ID columns", {
	mesh = make_mesh(.test_cells)
	expect_s3_class(mesh$pts, "data.table")
	expect_true(all(c("X", "Y", "ORIG_ID") %in% names(mesh$pts)))
})

test_that("make_mesh pruning reduces cell count", {
	mesh = make_mesh(.test_cells)
	expect_lte(nrow(mesh$pts), .test_N)
	expect_gt(nrow(mesh$pts), 0L)
})

test_that("make_mesh ORIG_ID values are valid indices into original cells", {
	mesh = make_mesh(.test_cells)
	expect_true(all(mesh$pts$ORIG_ID >= 1L))
	expect_true(all(mesh$pts$ORIG_ID <= .test_N))
	expect_equal(length(unique(mesh$pts$ORIG_ID)), nrow(mesh$pts))  # no duplicates
})

test_that("make_mesh edges is data.table with expected columns", {
	mesh = make_mesh(.test_cells)
	expect_s3_class(mesh$edges, "data.table")
	expected_cols = c("from_pt", "to_pt", "from_tri", "to_tri",
	                  "x0_pt", "x1_pt", "y0_pt", "y1_pt",
	                  "x0_tri", "x1_tri", "y0_tri", "y1_tri",
	                  "length_pt", "length_tri")
	expect_true(all(expected_cols %in% names(mesh$edges)))
})

test_that("make_mesh triangles is data.table with external column", {
	mesh = make_mesh(.test_cells)
	expect_s3_class(mesh$triangles, "data.table")
	expect_true("external" %in% names(mesh$triangles))
	expect_true(any(mesh$triangles$external))    # boundary tris were added
	expect_true(any(!mesh$triangles$external))   # real tris remain
})

test_that("make_mesh tri_to_pt is dgCMatrix with dims (num_tris x num_pts)", {
	mesh = make_mesh(.test_cells)
	expect_s4_class(mesh$tri_to_pt, "dgCMatrix")
	expect_equal(nrow(mesh$tri_to_pt), nrow(mesh$triangles))
	expect_equal(ncol(mesh$tri_to_pt), nrow(mesh$pts))
})

test_that("make_mesh adj is square dgCMatrix with dim num_pts", {
	mesh = make_mesh(.test_cells)
	expect_s4_class(mesh$adj, "dgCMatrix")
	expect_equal(nrow(mesh$adj), nrow(mesh$pts))
	expect_equal(ncol(mesh$adj), nrow(mesh$pts))
})

test_that("make_mesh adj is symmetric", {
	mesh = make_mesh(.test_cells)
	diff = mesh$adj - Matrix::t(mesh$adj)
	expect_equal(Matrix::nnzero(diff), 0L)
})

test_that("make_mesh weights is numeric vector of length num_edges", {
	mesh = make_mesh(.test_cells)
	expect_type(mesh$weights, "double")
	expect_equal(length(mesh$weights), nrow(mesh$edges))
})

test_that("make_mesh morse is NULL", {
	mesh = make_mesh(.test_cells)
	expect_null(mesh$morse)
})

test_that("make_mesh meta is an empty list", {
	mesh = make_mesh(.test_cells)
	expect_type(mesh$meta, "list")
	expect_length(mesh$meta, 0L)
})

test_that("make_mesh works when cells has no counts", {
	cells_no_counts = make_cells(.test_coords, .test_embeddings)
	mesh = make_mesh(cells_no_counts)
	expect_s3_class(mesh$pts, "data.table")
	expect_lte(nrow(mesh$pts), .test_N)
})

# ── Numerical equivalence tests ───────────────────────────────────────────────

test_that("make_mesh pts match old pipeline pts (X, Y, ORIG_ID)", {
	mesh = make_mesh(.test_cells)
	expect_equal(nrow(mesh$pts),      nrow(.test_mesh_old$pts))
	expect_equal(mesh$pts$X,          .test_mesh_old$pts$X)
	expect_equal(mesh$pts$Y,          .test_mesh_old$pts$Y)
	expect_equal(mesh$pts$ORIG_ID,    .test_mesh_old$pts$ORIG_ID)
})

test_that("make_mesh edges match old pipeline edges (from_pt, to_pt)", {
	mesh = make_mesh(.test_cells)
	expect_equal(nrow(mesh$edges),        nrow(.test_mesh_old$edges))
	expect_equal(mesh$edges$from_pt,      .test_mesh_old$edges$from_pt)
	expect_equal(mesh$edges$to_pt,        .test_mesh_old$edges$to_pt)
	expect_equal(mesh$edges$length_pt,    .test_mesh_old$edges$length_pt)
})

test_that("make_mesh triangles match old pipeline tris (X, Y, area, external)", {
	mesh = make_mesh(.test_cells)
	expect_equal(nrow(mesh$triangles),        nrow(.test_mesh_old$tris))
	expect_equal(mesh$triangles$X,            .test_mesh_old$tris$X)
	expect_equal(mesh$triangles$Y,            .test_mesh_old$tris$Y)
	expect_equal(mesh$triangles$area,         .test_mesh_old$tris$area)
	expect_equal(mesh$triangles$external,     .test_mesh_old$tris$external)
})

test_that("make_mesh tri_to_pt matches old pipeline tri_to_pt", {
	mesh = make_mesh(.test_cells)
	expect_equal(mesh$tri_to_pt, .test_mesh_old$tri_to_pt)
})

test_that("make_mesh adj connectivity matches old pipeline edge list", {
	mesh = make_mesh(.test_cells)
	# Every edge (from_pt, to_pt) must be present in adj
	edge_pairs = as.matrix(.test_mesh_old$edges[, .(from_pt, to_pt)])
	for (k in seq_len(min(nrow(edge_pairs), 50L))) {  # spot-check first 50
		i = edge_pairs[k, 1]
		j = edge_pairs[k, 2]
		expect_gt(mesh$adj[i, j], 0)
	}
})
