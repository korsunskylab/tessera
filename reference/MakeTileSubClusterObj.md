# Create a Seurat object for sub-clustering a specific cluster from matching tile and cell data

Create a Seurat object for sub-clustering a specific cluster from
matching tile and cell data

## Usage

``` r
MakeTileSubClusterObj(
  tile_obj,
  cluster,
  cell_obj,
  clusters.name = "seurat_clusters",
  npcs = 30,
  tile_assay = NULL,
  cell_assay = NULL,
  fast_sgd = TRUE,
  n_neighbors = 15,
  scale.factor = NULL,
  use.existing.embeddings = NULL,
  smooth_emb = c(0, 1),
  graph.name.cells = "cell_adj",
  meta.vars.include = NULL,
  harmony.group.by.vars = NULL,
  early_stop = TRUE
)
```

## Arguments

- tile_obj:

  A Seurat object containing tile-level data.

- cluster:

  The cluster label to sub-cluster.

- cell_obj:

  A Seurat object containing cell-level data with a `tile_id` metadata
  column.

- clusters.name:

  Name of the metadata column in `tile_obj` containing the cluster
  labels.

- npcs:

  Number of principal components to use for sub-clustering.

- tile_assay:

  Assay in `tile_obj` to use for sub-clustering. If NULL, use
  DefaultAssay.

- cell_assay:

  Assay in `cell_obj` to use for computing cell embeddings. If NULL, use
  DefaultAssay.

- fast_sgd:

  Whether to use fast SGD in UMAP.

- n_neighbors:

  Number of neighbors to use in UMAP.

- scale.factor:

  Scale factor for normalization. If NULL, use median nCount_RNA.

- use.existing.embeddings:

  Name of existing dimensional reduction in `tile_obj` to use for
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

## Value

A Seurat object for the tile subset that can be used for sub-clustering.
