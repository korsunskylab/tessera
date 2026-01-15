# Applies Tessera on a Seurat object

Applies Tessera on a Seurat object

## Usage

``` r
# S3 method for class 'Seurat'
GetTiles(
  obj,
  spatial,
  embeddings = NULL,
  dims.use = NULL,
  assay = NULL,
  group.by = NULL,
  raw_results = FALSE,
  tile.id.name = "tile_id",
  reduction.name = "pca",
  graph.name = "tile_adj",
  add.isolated.cells = TRUE,
  ...
)
```

## Arguments

- obj:

  Seurat object with spatial coordinates (and optionally, pre-computed
  single cell embeddings) stored as dimensional reductions.

- spatial:

  Name of dimensional reduction where the cells' x/y coordinates are
  stored.

- embeddings:

  Name of dimensional reduction where pre-computed single-cell
  embeddings are stored (a `num_cells` x `num_dim` matrix of cell
  embeddings across all latent dimensions). If missing, cell embeddings
  are calculated using PCA. If provided, the `npcs` parameter is
  ignored.

- assay:

  Seurat assay to pull data for when using the cell counts. Defaults to
  the DefaultAssay.

- group.by:

  Name of column in `obj@meta.data` to use for grouping cells into
  separate samples.

- raw_results:

  Whether to return the raw results from
  [`GetTiles.default()`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md).

- tile.id.name:

  Name of variable to store the tile IDs in the cell-level Seurat
  object.

- reduction.name:

  Name of dimensional reduction to store the aggregated tile-level
  embeddings in the tile-level Seurat object.

- graph.name:

  Name of graph to store tile adjacency matrix in the tile-level Seurat
  object.

- add.isolated.cells:

  Whether to add back isolated single cells that were pruned out. Only
  applies when `embeddings` are provided. Defaults to TRUE.

- ...:

  Arguments passed on to
  [`GetTiles.default`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)

  `meta_vars_include`

  :   Names of columns in meta_data to include in `dmt$pts`.

  `npcs`

  :   Number of PCs to compute for input to segmentation. Ignored if
      `embeddings` are provided directly.

  `smooth_emb`

  :   Number of smoothing iterations to perform on the cell embeddings
      prior to gradient computation. If a vector, then embeddings after
      each specified iteration are concatenated. If `0` is included,
      then the original embeddings are also included.

  `prune_thresh_quantile`

  :   Floating point value between 0 and 1, inclusive. Quantile of edge
      length above which edges are pruned. Defaults to 0.95.

  `prune_min_cells`

  :   Minimum number of cells required for a connected component of
      triangles to be kept. Defaults to 10.

  `prune_thresh`

  :   Edge length above which edges are pruned. If equal to NA, then
      this value is ignored and `thresh_quantile` is used to compute the
      threshold. Otherwise, if `thresh` is set, then `thresh_quantile`
      is ignored. Defaults to NA.

  `smooth_distance`

  :   One of `c('none', 'euclidean', 'projected', 'constant')`. If
      either `smooth_distance` or `smooth_similarity` is `'none'`, then
      no smoothing of the gradient field is conducted. Defaults to
      `'projected'`.

  `smooth_similarity`

  :   One of `c('none', 'euclidean', 'projected', 'constant')`. If
      either `smooth_distance` or `smooth_similarity` is `'none'`, then
      no smoothing of the gradient field is conducted. Defaults to
      `'projected'`.

  `smooth_iter`

  :   Number of rounds of gradient smoothing.

  `max_npts`

  :   Maximum number of cells allowed in each tile during the
      agglomerative clustering phase.

  `min_npts`

  :   Minimum number of cells allowed in each tile during the
      agglomerative clustering phase.

  `alpha`

  :   Parameter for scoring transcriptional similarity between adjacent
      tiles during the agglomerative clustering phase. For `alpha`, 0.2
      = conservative merging, 2 = liberal merging.

  `future.globals.maxSize`

  :   Maximum allowed size (in bytes) of global variables that are
      exported to each parallel worker. Increase this value if you get
      an error about global object size. Default is 8\*1024^3 (8 GB).

  `consolidate`

  :   Whether to consolidate results from multiple groups into a single
      collection of points and tiles (TRUE) or to return a list of
      separate results for each group (FALSE).

  `verbose`

  :   Whether to print progress messages for each stage of the
      segmentation pipeline.

## Value

A List containing a pair of Seurat objects:

1.  `obj`: the input single-cell object whose meta.data has been updated
    with tile assignments for each cell

2.  `tile_obj`: a Seurat object where each item represents an individual
    Tessera tile

## See also

Other GetTiles:
[`GetTiles()`](https://korsunskylab.github.io/tessera/reference/GetTiles.md),
[`GetTiles.default()`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)
