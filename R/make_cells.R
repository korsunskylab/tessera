#' Assemble cell-level data into a cells object
#'
#' Validates and structures raw inputs into the canonical cells list used
#' throughout the tessera pipeline. Does not compute embeddings — call
#' [do_pca()] or supply pre-computed embeddings before calling this function.
#'
#' @param coords N x 2 numeric matrix of cell coordinates. Columns are
#'   renamed to `x` and `y`.
#' @param embeddings Named list of N x D numeric matrices, one per embedding
#'   (e.g. `list(pca = do_pca(counts, npcs=10)$embeddings)`).
#' @param counts G x N [dgCMatrix][Matrix::dgCMatrix-class] of raw transcript
#'   counts. Optional; coerced to dgCMatrix if a different sparse format is
#'   supplied.
#' @param meta [data.table][data.table::data.table] with N rows of per-cell
#'   metadata. Optional; coerced from data.frame if needed.
#' @param params Named list of parameters carried alongside the object.
#'
#' @returns A cells list with slots:
#' * `coords`     — N x 2 numeric matrix
#' * `embeddings` — named list of N x D matrices
#' * `counts`     — G x N dgCMatrix (or NULL)
#' * `meta`       — data.table with N rows (or NULL)
#' * `params`     — list
#'
#' @export
make_cells = function(coords, embeddings, counts = NULL, meta = NULL, params = list()) {
	if (!is.matrix(coords) || ncol(coords) != 2)
		stop("coords must be an N x 2 numeric matrix")
	N = nrow(coords)
	colnames(coords) = c("x", "y")

	if (!is.null(embeddings)) {
		if (!is.list(embeddings) || is.null(names(embeddings)))
			stop("embeddings must be a named list")
		for (nm in names(embeddings)) {
			if (!is.matrix(embeddings[[nm]]) || nrow(embeddings[[nm]]) != N)
				stop(sprintf("embeddings$%s must be a matrix with N rows", nm))
		}
	}

	if (!is.null(counts)) {
		if (!inherits(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
		if (ncol(counts) != N)
			stop("counts must have N columns (one per cell)")
	}

	if (!is.null(meta)) {
		if (!is.data.table(meta)) meta = as.data.table(meta)
		if (nrow(meta) != N)
			stop("meta must have N rows (one per cell)")
	}

	list(
		coords     = coords,
		embeddings = embeddings,
		counts     = counts,
		meta       = meta,
		params     = params
	)
}
