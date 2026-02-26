# Shared test fixtures for tessera unit tests.
# Loaded automatically by testthat before any test file runs.
# Each block adds its own golden reference outputs here as the refactor
# progresses. Old-API calls are replaced block by block.

library(data.table)
data(tessera_warmup)

# ── Raw inputs ────────────────────────────────────────────────────────────────
.test_X      = tessera_warmup$meta_data$X
.test_Y      = tessera_warmup$meta_data$Y
.test_counts = tessera_warmup$counts          # 479 genes x 3177 cells dgCMatrix
.test_meta   = as.data.table(tessera_warmup$meta_data)
.test_N      = length(.test_X)               # 3177

# coords as N x 2 matrix (new spec format)
.test_coords = matrix(
	c(.test_X, .test_Y),
	ncol = 2,
	dimnames = list(NULL, c("x", "y"))
)

# ── Block 1 golden reference: PCA embeddings ──────────────────────────────────
.test_pca        = do_pca(.test_counts, npcs = 10)
.test_embeddings = list(pca = .test_pca$embeddings)

# Pre-built cells object reused across blocks (avoids re-running PCA each time)
.test_cells = make_cells(.test_coords, .test_embeddings, .test_counts, .test_meta)

# ── Block 2 golden reference: mesh from old pipeline ──────────────────────────
.test_mesh_old = local({
	d = init_data(.test_X, .test_Y, .test_counts)
	d = prune_graph(d)
	add_exterior_triangles(d)
})

# ── Block 3 golden reference: field from old pipeline ─────────────────────────
# PCA is run on PRUNED cells only, matching the old GetTiles.default pipeline.
.test_mesh        = make_mesh(.test_cells)
.test_pca_pruned  = do_pca(.test_counts[, .test_mesh$pts$ORIG_ID], npcs = 10)

# cells object with embeddings aligned to the pruned mesh (nrow = nrow(mesh$pts))
.test_cells_field = local({
	cells            = .test_cells
	cells$embeddings = list(pca = .test_pca_pruned$embeddings)
	cells
})

# old-pipeline field output (compute_gradients + compress_gradients_svd)
.test_field_old = local({
	dmt            = .test_mesh_old
	dmt$udv_cells  = list(embeddings = .test_pca_pruned$embeddings)
	field = compute_gradients(dmt,
		smooth_distance  = 'projected',
		smooth_similarity = 'projected',
		smooth_iter      = 1)
	compress_gradients_svd(field)
})

# ── Block 4 golden reference: morse from old pipeline ─────────────────────────
# (also used as the base for Block 5)
# data.table::copy prevents dmt_set_f from mutating .test_mesh_old slots
.test_morse_old = local({
	dmt = list(
		pts   = data.table::copy(.test_mesh_old$pts),
		tris  = data.table::copy(.test_mesh_old$tris),
		edges = data.table::copy(.test_mesh_old$edges)
	)
	dmt       = dmt_set_f(dmt, .test_field_old)
	dmt$prim  = do_primary_forest(dmt)
	dmt$dual  = do_dual_forest(dmt)
	dmt$e_sep = dmt_get_separatrices(dmt)
	dmt
})

# ── Block 5 golden reference: tile assignment from old separatrices ───────────
# dmt_assign_tiles uses old igraph API (broken in igraph 2.x); build the
# equivalent reference with the modern make_empty_graph + add_edges approach.
.test_dmt_old = local({
	n_pts    = nrow(.test_morse_old$pts)
	n_edges  = nrow(.test_morse_old$edges)
	e_keep   = setdiff(seq_len(n_edges), .test_morse_old$e_sep)
	edge_vec = c(rbind(
		.test_morse_old$edges$from_pt[e_keep],
		.test_morse_old$edges$to_pt[e_keep]
	))
	g = igraph::make_empty_graph(n = n_pts, directed = FALSE)
	if (length(edge_vec) > 0) g = igraph::add_edges(g, edge_vec)
	list(agg_id = igraph::components(g)$membership)
})

# mesh with morse populated — reused in Block 5 tests to avoid recomputing
.test_mesh_morse = local({
	m        = .test_mesh
	m$morse  = compute_morse(compute_field(.test_cells_field, .test_mesh), .test_mesh)
	m
})
