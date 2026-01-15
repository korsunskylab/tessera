# Generic function that runs the Tessera algorithm on single-cell spatial data

GetTiles is a generic function that runs the main Tessera algorithm. If
working with a Seurat object, please refer to the documentation of the
appropriate generic API:
[`GetTiles.Seurat()`](https://korsunskylab.github.io/tessera/reference/GetTiles.Seurat.md).
If users work with other forms of the input, they can pass them directly
to Tessera using the
[`GetTiles.default()`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)
API. The function arguments listed here are common in all GetTiles
interfaces.

## Usage

``` r
GetTiles(...)
```

## Arguments

- ...:

  Arguments passed on to
  [`GetTiles.default`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)

  `meta_vars_include`

  :   Names of columns in meta_data to include in `dmt$pts`.

  `group.by`

  :   Name of column in `meta_data` that provides the group IDs. Tessera
      tiles are constructed separately for each group (which could be
      separate experimental samples or FOVs).

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

If used with a Seurat object, it will return a pair of Seurat objects:

1.  the input single-cell object updated with tile assignments for each
    cell, and

2.  a Seurat object where each item represents an individual Tessera
    tile.

For standalone operation, it returns Lists with the output of Tessera
segmentation (see
[`GetTiles.default()`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)).

## See also

Other GetTiles:
[`GetTiles.Seurat()`](https://korsunskylab.github.io/tessera/reference/GetTiles.Seurat.md),
[`GetTiles.default()`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)
