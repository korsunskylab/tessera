# Updates information for boundaries after merging two tiles

For every boundary that exists between the new merged tile and other
tiles, the following values must be updated:

- `E_pcs_merge[e_update,]`: PCs that would result from future merging.

- `E_perimeter_merge[e_update]`: Perimeters that would result from
  future merging.

- `E_w[e_update]`: Expression similarity score.

- `E_score_size[e_update]`: Size of merged tile score.

- `dC[e_update]`: Delta shape compactness score for merged tile.

- `E_dscore[e_update]`: Overall merging score (product of `w`,
  `score_size`, `dC`) Note that `E_npts` and `E_area` should already
  have been previously updated.

## Usage

``` r
update_E_cpp(
  V_pcs,
  V_perimeter,
  V_area,
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
  e_update,
  V_to_E_from,
  V_to_E_to,
  d_mu,
  d_sig,
  agg_mode,
  min_npts,
  max_npts
)
```

## Arguments

- e_update:

  Edges that should be updated.
