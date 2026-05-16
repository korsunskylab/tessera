# Create a Seurat object for sub-clustering a specific cluster

Create a Seurat object for sub-clustering a specific cluster

## Usage

``` r
MakeSubClusterObj(
  obj,
  cluster,
  clusters.name = "seurat_clusters",
  npcs = 30,
  fast_sgd = TRUE,
  n_neighbors = 15,
  assay = NULL,
  scale.factor = NULL,
  use.existing.embeddings = NULL,
  meta.vars.include = NULL,
  harmony.group.by.vars = NULL,
  early_stop = TRUE,
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

- npcs:

  Number of principal components to use for sub-clustering.

- fast_sgd:

  Whether to use fast SGD in UMAP.

- n_neighbors:

  Number of neighbors to use in UMAP.

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

- ...:

  Additional parameters to pass to
  [`harmony::RunHarmony`](https://pati-ni.github.io/harmony/reference/RunHarmony.html).
  Ignored if `harmony.group.by.vars` is NULL.

## Value

A Seurat object for the sub-clustered cells.
