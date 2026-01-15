# Merges tiles using single-linkage agglomerative clustering

Note that this function mutates many of the input values to keep track
of updated values for tiles and their borders after each successive
merge. Every merge of adjacent tiles has an associated score (higher
means merging is more favorable). Merging is conducted greedily, one
step of a time, updating the score associated with pair of tiles at each
step.

## Usage

``` r
merge_aggs(
  aggs,
  agg_mode,
  d_mu = NULL,
  d_sig = NULL,
  iter_max = NULL,
  dscore_thresh = 0,
  min_npts = 0,
  max_npts = Inf,
  min_area = 0,
  max_area = Inf
)
```

## Value

The `aggs` data structure with updated values for

- meta_data:

  A data table with attributes for the new merged tiles: \* `X,Y`:
  Centroid of each tile. (NOT updated) \* `npts`: Number of points in
  each tile. (Updated) \* `shape`: A `sfc` list-column with the
  geometries for each tile. (NOT updated) \* `area,perimeter`: Area and
  perimeter of each tile. (Updated)

- edges:

  A data table with attributes for edges between adjacent merged tiles,
  where edges with infinite `dscore` are removed: \* `from,to`: Tile IDs
  for the two tiles bordering this edge. (Updated for new merged tiles)
  \* `x0,y0,x1,y1`: Centroid coordinates for the two tiles bordering
  this edge. (NOT updated) \* `area,npts`: Sum of areas and numbers of
  points in the two tiles bordering this edge. (Updated) \*
  `edge_length`: Total length of the border between the `from` and `to`
  tiles. (Updated) \* `dscore`: Overall score for merging two tiles.
  Product of `w`, `score_size`, and `dC`. (Updated) \* `w`: Gene
  expression similarity score. (Updated) \* `score_size`: Penalizes
  tiles with many points. (Updated) \* `perimeter_merge`: Perimeter of
  merged tile. (Updated)

- pcs:

  A `num_tiles` x `npcs` matrix with the average embedding value over
  all cells in each tile. (Updated)

- pcs_merged:

  A `num_edges` x `npcs` matrix with average PCs for merged tiles.
  (Updated)

- d_mu,d_sig:

  Parameters used to calculate `w`.

- aggmap:

  A length `orig_num_tiles` vector mapping each original tile ID to the
  new tile IDs after merging.
