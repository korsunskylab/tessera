# Find sub-clusters within a specific cluster of a Seurat object using matching tile and cell data

Find sub-clusters within a specific cluster of a Seurat object using
matching tile and cell data

## Usage

``` r
FindTileSubCluster(
  obj,
  cluster,
  cell_obj = NULL,
  clusters.name = "seurat_clusters",
  tile_assay = NULL,
  cell_assay = NULL,
  sub.clusters.name = NULL,
  resolution = 0.2,
  algorithm = 4,
  npcs = 20,
  method = "igraph",
  n_neighbors = 25,
  fast_sgd = TRUE,
  scale.factor = NULL,
  use.existing.embeddings = NULL,
  smooth_emb = c(0, 1),
  graph.name.cells = "cell_adj",
  meta.vars.include = NULL,
  harmony.group.by.vars = NULL,
  early_stop = TRUE,
  return_obj_sub = FALSE
)
```

## Arguments

- obj:

  A Seurat object containing tile-level data.

- cluster:

  The cluster label to sub-cluster.

- cell_obj:

  (Opional) A Seurat object containing cell-level data with a `tile_id`
  metadata column. If NULL, sub-clustering is performed using only `obj`
  without recomputing tile embeddings from cell data.

- clusters.name:

  Name of the metadata column in `obj` containing the cluster labels.

- tile_assay:

  Assay in `obj` to use for sub-clustering. If NULL, use DefaultAssay.

- cell_assay:

  Assay in `cell_obj` to use for computing cell embeddings. If NULL, use
  DefaultAssay.

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
  sub-clustering. If NULL, compute PCA on the subsetted cell-level data.

- smooth_emb:

  Number of smoothing iterations to perform on the cell embeddings, as
  for `GetTiles`. If a vector, then embeddings after each specified
  iteration are concatenated. If `0` is included, then the original
  embeddings are also included.

- graph.name.cells:

  Name of the graph in `cell_obj` containing the cell adjacency graph.

- meta.vars.include:

  Metadata variables to include in the sub-clustered object.

- harmony.group.by.vars:

  Metadata variables to use for Harmony integration.

- early_stop:

  Whether to use early stopping in Harmony.

- return_obj_sub:

  If TRUE, return a list with the updated `obj` and the sub-clustered
  object.

## Value

The input Seurat object with sub-cluster labels added to metadata.
