# Add sub-cluster labels to a Seurat object

Add sub-cluster labels to a Seurat object

## Usage

``` r
AddSubClusterLabels(
  obj,
  obj_sub,
  cluster,
  clusters.name = "seurat_clusters",
  sub.clusters.name = NULL
)
```

## Arguments

- obj:

  A Seurat object.

- obj_sub:

  A Seurat object containing the sub-clustered cells.

- cluster:

  The parent cluster label in `obj` that was sub-clustered.

- clusters.name:

  Name of the metadata column in `obj` containing the parent cluster
  labels.

- sub.clusters.name:

  Name of the metadata column to store the sub-cluster labels. If NULL,
  defaults to `<clusters.name>_s<cluster>`.

## Value

The input Seurat object with sub-cluster labels added to metadata.
