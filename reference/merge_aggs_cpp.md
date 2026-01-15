# Merges tiles using single-linkage agglomerative clustering

Note that this function mutates many of the inputs to keep track of
updated values for tiles and their borders after each successive merge.
Every merge of adjacent tiles has an associated score (higher means
merging is more favorable). Merging is conducted greedily, one step of a
time, updating the score associated with pair of tiles at each step.

## Usage

``` r
merge_aggs_cpp(
  V_pcs,
  V_area,
  V_perimeter,
  V_npts,
  E_from,
  E_to,
  E_npts,
  E_area,
  E_edge_length,
  E_pcs_merge,
  E_w,
  E_perimeter_merge,
  E_score_size,
  E_dscore,
  d_mu,
  d_sig,
  iter_max,
  agg_mode,
  dscore_thresh,
  min_npts,
  max_npts
)
```

## Arguments

- V_pcs, V_area, V_perimeter, V_npts:

  Metadata associated with each each tile. (Updated)

- E_from, E_to:

  Length `num_edges` vectors associated with each pair of adjacent
  tiles. Specifies the two tiles that each border (0-indexed). (Updated)

- E_npts, E_area, E_edge_length, E_pcs_merge:

  Metadata associated with each pair of adjacent tiles. (Updated)

- E_w, E_perimeter_merge, E_score_size, E_dscore:

  Scores associated with merging each pair of adjacent tiles. (Updated)

## Value

A list of `orig_num_tiles` vectors, disjoint sets specifying the IDs for
which of the original tiles have been merged together. Some vectors will
have length 0.
