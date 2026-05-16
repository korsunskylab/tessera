# Find sub-clusters within a specific cluster of a Seurat object

Find sub-clusters within a specific cluster of a Seurat object

## Usage

``` r
FindSubClusterCustom(
  obj,
  cluster,
  clusters.name = "seurat_clusters",
  sub.clusters.name = NULL,
  resolution = 0.5,
  algorithm = 1,
  npcs = 30,
  method = "igraph",
  n_neighbors = 15,
  fast_sgd = TRUE,
  scale.factor = NULL,
  use.existing.embeddings = NULL,
  meta.vars.include = NULL,
  harmony.group.by.vars = NULL,
  early_stop = TRUE,
  return_obj_sub = FALSE,
  ...
)
```

## Arguments

- obj:

  A Seurat object.

- cluster:

  The cluster label to sub-cluster.

- clusters.name:

  Name of the metadata column in `obj` containing the cluster labels.

- sub.clusters.name:

  Name of the metadata column to store the sub-cluster labels. If NULL,
  defaults to `<clusters.name>_s<cluster>`.

- resolution:

  Resolution parameter for clustering.

- algorithm:

  Clustering algorithm to use. See
  [`?Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
  for details.

- npcs:

  Number of principal components to use for sub-clustering.

- method:

  Method to use for clustering. See
  [`?Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)
  for details.

- n_neighbors:

  Number of neighbors to use in UMAP.

- fast_sgd:

  Whether to use fast SGD in UMAP.

- scale.factor:

  Scale factor for normalization. If NULL, use median nCount_RNA.

- use.existing.embeddings:

  Name of existing dimensional reduction in `obj` to use for
  sub-clustering. If NULL, compute PCA on the subsetted data.

- meta.vars.include:

  Metadata variables to include in the sub-clustered object.

- harmony.group.by.vars:

  Metadata variables to use for Harmony integration.

- early_stop:

  Whether to use early stopping in Harmony.

- return_obj_sub:

  If TRUE, return a list with the updated `obj` and the sub-clustered
  object.

- ...:

  Additional parameters to pass to
  [`harmony::RunHarmony`](https://pati-ni.github.io/harmony/reference/RunHarmony.html).
  Ignored if `harmony.group.by.vars` is NULL.

## Value

The input Seurat object with sub-cluster labels added to metadata.
