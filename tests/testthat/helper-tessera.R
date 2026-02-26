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
