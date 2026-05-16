# Merge clusters in a Seurat object

Merge clusters in a Seurat object

## Usage

``` r
MergeClusters(
  obj,
  to_merge,
  new_label,
  clusters.name = "seurat_clusters",
  new.clusters.name = NULL
)
```

## Arguments

- obj:

  A Seurat object.

- to_merge:

  A vector of cluster labels to merge.

- new_label:

  The new label for the merged cluster.

- clusters.name:

  Name of the metadata column in `obj` containing the cluster labels to
  be merged. Defaults to 'seurat_clusters'.

- new.clusters.name:

  Name of the metadata column to store the new cluster labels. If NULL,
  defaults to `<clusters.name>_m<new_label>`.

## Value

The input Seurat object with the merged cluster labels added to
metadata.
