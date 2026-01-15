# Run UMAP and save fgraph and embeddings in Seurat object

Run UMAP and save fgraph and embeddings in Seurat object

## Usage

``` r
RunUMAPCustom(
  obj,
  reduction = "pca",
  dims = NULL,
  fgraph_only = FALSE,
  graph.name = NULL,
  reduction.name = "umap",
  assay = NULL,
  key = "UMAP_",
  n_neighbors = 30,
  n_components = 2,
  metric = "cosine",
  spread = 1,
  min_dist = 0.3,
  n_threads = NULL,
  fast_sgd = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- obj:

  A Seurat object.

- reduction:

  Name of dimensional reduction in `obj` to use as input to UMAP.

- dims:

  Dimensions of `reduction` to use as input to UMAP. If NULL, use all
  dimensions.

- fgraph_only:

  If TRUE, only compute and store the fuzzy simplicial graph (fgraph)
  and skip UMAP embedding.

- graph.name:

  Name of graph to store the fgraph in `obj`. Defaults to
  `<assay>_fgraph`.

- reduction.name:

  Name of dimensional reduction to store the UMAP embeddings in `obj`.
  Defaults to "umap".

- assay:

  Assay to set as default assay for the new dimensional reduction. If
  NULL, use the default assay of `obj`.

- key:

  Key prefix to use for the new dimensional reduction. Defaults to
  "UMAP\_".

- n_neighbors:

  Number of nearest neighbors to use in UMAP. See
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for details.

- n_components:

  Number of UMAP dimensions to compute. Ignored if `fgraph_only` is
  TRUE.

- metric:

  Distance metric to use in UMAP. See
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for details.

- spread:

  UMAP spread parameter. See
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for details.

- min_dist:

  UMAP minimum distance parameter. See
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for details.

- n_threads:

  Number of threads to use in UMAP. See
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for details.

- fast_sgd:

  Whether to use the fast stochastic gradient descent optimization in
  UMAP. See
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for details.

- verbose:

  Whether to print progress messages.

- ...:

  Additional parameters to pass to
  [`uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

The input Seurat object with the fgraph and (optionally) UMAP embeddings
added.
