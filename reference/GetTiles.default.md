# Run full DMT segmentation pipeline to make aggregated tiles from cells

Segmentation has four main steps:

1.  **Preparing data structures:** A triangle mesh is constructed using
    Delauney triangulation and pruned to eliminate long edges. PC
    embeddings for each each are computed.

2.  **Computing gradients:** Gradients are calculated at each point by
    considering the difference in expression between each cell in its
    neighbors in the mesh. These gradients are smoothed using
    (anisotropic) bilateral filtering, and then gradients are defined
    for edges and triangles in the mesh by averaging the points that
    each edge or triangle contains.

3.  **DMT:** A scalar field is defined by taking the magnitude of the
    total gradient at each point/edge/triangle. Then DMT-based
    segmentation is performed by constructing a maximum spanning forest
    on the triangles and a minimum spanning forest on the points.
    Separatrices that separate cells into tiles of homogeneous
    composition are defined by tracing paths between critical points,
    particularly between saddle edges and maximum triangles.

4.  **Aggregation:** Tiles from DMT-based segmentation are merged using
    single-linkage agglomerative clustering to obtain tiles containing a
    number of cells between a user-provided minimum and maximum value.
    Pairs of adjacent tiles are scored according to their
    transcriptional similarity, compactness of shape after merging, and
    number of cells in order to prioritize favorable merges in each
    agglomerative clustering step.

## Usage

``` r
# Default S3 method
GetTiles(
  X,
  Y,
  counts = NULL,
  embeddings = NULL,
  loadings = NULL,
  meta_data = NULL,
  meta_vars_include = NULL,
  group.by = NULL,
  npcs = 25,
  smooth_emb = c(0, 1),
  prune_thresh_quantile = 0.99,
  prune_min_cells = 1,
  prune_thresh = NA,
  smooth_distance = c("none", "euclidean", "projected", "constant")[3],
  smooth_similarity = c("none", "euclidean", "projected", "constant")[3],
  smooth_iter = 1,
  on_edges = FALSE,
  max_npts = 50,
  min_npts = 5,
  alpha = 1,
  .progress = TRUE,
  .options = NULL,
  future.globals.maxSize = 8 * 1024^3,
  consolidate = TRUE,
  verbose = FALSE
)
```

## Arguments

- X, Y:

  A pair of numeric vectors with the coordinates for each of `num_cells`
  points.

- counts:

  A `num_genes` x `num_cells` gene-by-cell matrix of transcript counts.
  Optional if `embeddings` are provided directly.

- embeddings:

  A `num_cells` x `num_dim` matrix of cell embeddings across all latent
  dimensions. If missing, cell embeddings are calculated using PCA. If
  provided, the `npcs` parameter is ignored.

- loadings:

  (Optional) A `num_genes` x `num_dim` matrix of gene loadings.

- meta_data:

  A data frame with additional cell metadata to include in `dmt$pts`.

- meta_vars_include:

  Names of columns in meta_data to include in `dmt$pts`.

- group.by:

  Name of column in `meta_data` that provides the group IDs. Tessera
  tiles are constructed separately for each group (which could be
  separate experimental samples or FOVs).

- npcs:

  Number of PCs to compute for input to segmentation. Ignored if
  `embeddings` are provided directly.

- smooth_emb:

  Number of smoothing iterations to perform on the cell embeddings prior
  to gradient computation. If a vector, then embeddings after each
  specified iteration are concatenated. If `0` is included, then the
  original embeddings are also included.

- prune_thresh_quantile:

  Floating point value between 0 and 1, inclusive. Quantile of edge
  length above which edges are pruned. Defaults to 0.95.

- prune_min_cells:

  Minimum number of cells required for a connected component of
  triangles to be kept. Defaults to 10.

- prune_thresh:

  Edge length above which edges are pruned. If equal to NA, then this
  value is ignored and `thresh_quantile` is used to compute the
  threshold. Otherwise, if `thresh` is set, then `thresh_quantile` is
  ignored. Defaults to NA.

- smooth_distance:

  One of `c('none', 'euclidean', 'projected', 'constant')`. If either
  `smooth_distance` or `smooth_similarity` is `'none'`, then no
  smoothing of the gradient field is conducted. Defaults to
  `'projected'`.

- smooth_similarity:

  One of `c('none', 'euclidean', 'projected', 'constant')`. If either
  `smooth_distance` or `smooth_similarity` is `'none'`, then no
  smoothing of the gradient field is conducted. Defaults to
  `'projected'`.

- smooth_iter:

  Number of rounds of gradient smoothing.

- max_npts:

  Maximum number of cells allowed in each tile during the agglomerative
  clustering phase.

- min_npts:

  Minimum number of cells allowed in each tile during the agglomerative
  clustering phase.

- alpha:

  Parameter for scoring transcriptional similarity between adjacent
  tiles during the agglomerative clustering phase. For `alpha`, 0.2 =
  conservative merging, 2 = liberal merging.

- future.globals.maxSize:

  Maximum allowed size (in bytes) of global variables that are exported
  to each parallel worker. Increase this value if you get an error about
  global object size. Default is 8\*1024^3 (8 GB).

- consolidate:

  Whether to consolidate results from multiple groups into a single
  collection of points and tiles (TRUE) or to return a list of separate
  results for each group (FALSE).

- verbose:

  Whether to print progress messages for each stage of the segmentation
  pipeline.

## Value

If `consolidate==TRUE`, a List with the results of segmentation,
combined across groups (otherwise, if `consolidate==FALSE`, a named
List, which contains separate results for each group):

- dmt:

  Mesh data structures with input points/edges/triangles and the results
  from segmentation:

  - `pts`: A data table with `num_cells_pruned` rows containing cells in
    the mesh that remain after Delauney triangulation and pruning, with
    the following columns:

    - `X`,`Y`: Coordinates of each cell.

    - `ORIG_ID`: Index of each cell in the original inputs to
      [`GetTiles()`](https://korsunskylab.github.io/tessera/reference/GetTiles.md)
      (`X`, `Y`, `counts`, `meta_data`)

    - Columns from `meta_vars_include`.

    - `f`: Scalar value used for initial DMT-based segmentation,
      computed from the spatial gradient in expression at each point.

    - `agg_id`: Unique ID for the tile that each point belongs to in the
      final segmentation.

  - `tris`: A data table with `num_triangles` rows containing triangles
    in the mesh after Delauney triangulation and pruning, with the
    following columns:

    - `X`,`Y`: Coordinates of each triangle's centroid.

    - `area`: Area of triangle.

    - `height`: Largest height of each triangle.

    - `external`: Logical value that is `TRUE` if the triangle is a
      degenerate triangle that was added along a boundary edge to ensure
      that every edge is adjacent to two triangles. Degenerate triangles
      only have two vertices (the endpoints of the boundary edge).

    - `f`: Scalar value used for initial DMT-based segmentation,
      computed from the spatial gradient in expression for triangle.

  - `edges`: A data table with `num_edges` rows containing edges between
    adjacent points and triangles in the mesh after Delauney
    triangulation and pruning, with the following columns:

    - `from_pt`,`to_pt`: Indices of the two adjacent cells that are
      connected by an edge.

    - `from_tri`,`to_tri`: Indices of the two adjacent triangles that
      are connected by an edge.

    - `x0_pt`,`x1_pt`,`y0_pt`,`y1_pt`: Coordinates of the two adjacent
      cells that are connected by an edge.

    - `x0_tri`, `x1_tri`, `y0_tri`, `y1_tri`: Coordinates of the two
      adjacent triangles (centroids) that are connected by an edge.

    - `length_pt`, `length_tri`: Distance between adjacent cells or
      triangle centroids. Used for pruning step. (Warning - these are
      computed prior to pruning and are not updated after pruning and
      adding exterior triangles.)

    - `boundary`: Logical value that is `TRUE` if the edge is at the
      boundary and is adjacent to only a single internal triangle. For
      boundary edges, a degenerate external triangle is added along the
      boundary edge to ensure that every edge is adjacent to two
      triangles.

    - `f_prim`,`f_dual`: Scalar values used for initial DMT-based
      segmentation, computed from the spatial gradient in expression at
      each edge. Primal edges connect points and average the `f` values
      at the two adjacent points. Dual edges connect triangles and
      average the `f` values at the two adjacent triangles.

    - `agg_from`,`agg_to`: Unique ID for the tile that each adjacent
      point belongs to in the final segmentation.

  - `tri_to_pt`: A `num_triangles` x `num_cells_pruned` sparse matrix
    with value 1 at `(i,j)` if triangle `i` has point `j` as a vertex.
    Each internal triangle has 3 vertices, and each degenerate external
    triangle has 2 vertices.

  - `counts`: A `num_genes` x `num_cells_pruned` gene-by-cell matrix of
    transcript counts.

  - `udv_cells`: A List with the PC embeddings for each cell, which are
    used for segmentation.

    - `loadings`: If not provided as input, a `num_genes` x `npcs`
      matrix of gene loadings for each PC. Each column is a unit vector.

    - `embeddings`: If not provided as input, a `num_cells` x `npcs`
      matrix of cell embeddings across all PCs. Each column `j` has
      magnitude equal to the `j`th singular value. That is, PCs with
      larger contribution to the total variance will have embeddings of
      proportionally larger magnitude.

  - `prim`: The primal minimum spanning forest on points. A List with
    the following attributes:

    - `edges`: A data table with `forest_size` rows, where each row is a
      directed edge in the minimum spanning forest. There are six
      columns:

      - `from,to`: Index of source and target points for each edge.

      - `x0,y0`: Coordinates of source point for each edge.

      - `x1,y1`: Coordinates of target point for each edge.

    - `saddles`: A length `num_saddles` vector with edge indices for
      possible saddle edges.

    - `labels`: A length `num_points` vector of labels for the connected
      components in the minimum spanning tree. Each connected component
      is labeled by the index of its critical point.

    - `minima`: A length `num_critpts` vector of critical points
      (minima).

    - `parent`: A length `num_points` vector containing the parent
      (source) point for each point in the directed spanning forest.
      Critical points have no parent, so the value should be ignored.

    - `parent_edge`: A length `num_points` vector containing the
      directed edge that has each point as a target node. Critical
      points have no parent edge, so the value should be ignored.

  - `dual`: The dual maximum spanning forest on triangles. A List with
    the following attributes:

    - `edges`: A data.table with `forest_size` rows, where each row is a
      directed edge in the maximum spanning forest. There are six
      columns:

      - `from,to`: Index of source and target triangles for each edge.

      - `x0,y0`: Coordinates of source triangle for each edge.

      - `x1,y1`: Coordinates of target triangle for each edge.

    - `saddles`: A length `num_saddles` vector with edge indices for
      possible saddle edges.

    - `labels`: A length `num_triangles` vector of labels for the
      connected components in the maximum spanning tree. Each connected
      component is labeled by the index of its critical triangle.

    - `maxima`: A length `num_critpts` vector of critical triangles
      (maxima).

    - `parent`: A length `num_triangles` vector containing the parent
      (source) triangle for each triangle in the directed spanning
      forest. Critical triangles have no parent, so the value should be
      ignored.

    - `parent_edge`: A length `num_triangles` vector containing the
      directed edge that has each triangle as a target node. Critical
      triangles have no parent edge, so the value should be ignored.

  - `e_sep`: A length `num_sep_edges` vector of edge indices that make
    up the separatrices, which separate points into different
    components.

- aggs:

  The tiles that result from DMT-based segmentation and agglomeration. A
  List data structure stores the tiles and their adjacencies using the
  following attributes:

  - `meta_data`: A data table with `num_tiles` rows with metadata for
    each tile:

    - `ID`: Unique ID for each tile.

    - `X`,`Y`: Centroid of each tile.

    - `npts`: Number of points in each tile.

    - `shape`: A `sfc` list-column with the geometries for each tile.

    - `area,perimeter`: Area and perimeter of each tile.

  - `edges`: Additional attributes are calculated:

    - `from,to`: Tile IDs for the two tiles bordering this edge.

    - `x0,y0,x1,y1`: Centroid coordinates for the two tiles bordering
      this edge.

    - `area,npts`: Sum of areas and numbers of points in the two tiles
      bordering this edge.

    - `edge_length`: Total length of the border between the `from` and
      `to` tiles.

    - `dscore`: Overall score for merging two tiles. Product of `w`,
      `score_size`, and `dC`.

    - `w`: Gene expression similarity score.

    - `score_size`: Penalizes tiles with many points.

    - `perimeter_merge`: Perimeter of merged tile.

  - `pcs`: A `num_tiles` x `npcs` matrix with the average embedding
    value over all cells in each tile.

  - `pcs_merged`: A `num_edges` x `npcs` matrix with average PCs for the
    new tile if the two adjacent tiles connected by the edge were
    merged.

  - `d_mu`,`d_sig`: Parameters used to calculate `w` in the edge score
    `dscore`.

  - `aggmap`: A length `orig_num_tiles` vector mapping each original
    tile ID to the new tile IDs after merging.

  - `adj`: Sparse adjacency matrix between all tiles (if
    `consolidate==TRUE`).

  - `counts`: A `num_genes` x `num_tiles` gene-by-tile matrix of
    aggregated transcript counts.

## Details

### Computing gradients (smoothing)

Gradient fields are smoothed using bilateral filtering, in which the
smoothed gradient of each point is computed as the weighted average of
the neighbors' gradients, considering both distance in space and also
similarity in gradients. The weight of each neighbor is computed from
the product of two scores:

- `distance` score: Generally, closer neighbors have greater weight.

  - if `'euclidean'`: Gaussian transformation of the Euclidean distance
    of a cell from its neighbor, so that more distant neighbors have
    less weight.

  - if `'projected'`: An anisotropic filter that accounts for expected
    change in expression along the direction of the neighbor. The
    expected change in expression is calculated from the gradient field
    as the total derivative in the direction of the neighbor. This
    change in expression is then Gaussian transformed so that neighbors
    that are more distant along the direction of greatest change have
    less weight.

  - if `'constant'`: All neighbors have equal `distance` weights

- `similarity` score: Generally, neighbors with more similar gradients
  have greater weight

  - if `'euclidean'`: Gaussian transformation of the Euclidean distance
    between a cell's gradient field and its neighbor's gradient field.

  - if `'projected'`: Gaussian transformation of the cosine distance
    between a cell's gradient field and its neighbor's gradient field.

  - if `'constant'`: All neighbors have equal `similarity` weights

### Aggregation (scores)

Pairs of adjacent tiles are scored according to their transcriptional
similarity, compactness of shape after merging, and number of cells in
order to prioritize favorable merges in each agglomerative clustering
step. The `dscore` for each edge is computed as the product of the
following three factors, where higher scores favor merging of adjacent
tiles:

- `w`: A 2-cluster GMM is used to determine the mean `mu` and standard
  deviation `sig` of the distance between tiles that have similar gene
  expression (Euclidean distance `d` in PC space). Then we define
  `d_mu = mu + sig` and `d_sig = alpha * sig`, and calculate `w` as
  `w = 0.5 - 1 / (1 + exp(-(d - d_mu) / d_sig))`. Ranges from -0.5 to
  0.5. If adjacent tiles are very dissimilar (`d >> d_mu`), then `d` is
  large, and `w` is close to `-0.5`. If adjacent tiles are very similar
  (`d < d_mu`), then `d` is small, and `w` is positive.

- `score_size`: `(1 - npts_from/max_npts) * (1 - npts_to/max_npts)`.
  Ranges from 0 to 1. In the first round of aggregation, if merging the
  two tiles would have a total number of points ≥`max_npts`, then
  `score_size` is set to `-Inf`, which prevents merging. In the second
  round of aggregation, `dscore` is set to `-Inf` if both adjacent tiles
  have at least `min_npts` cells, to prioritize merging of small tiles.

- `dC`: `.5 * (C_merge - C_from - C_to + 1)`. Ranges from 0 to 1.

## See also

Other GetTiles:
[`GetTiles()`](https://korsunskylab.github.io/tessera/reference/GetTiles.md),
[`GetTiles.Seurat()`](https://korsunskylab.github.io/tessera/reference/GetTiles.Seurat.md)
