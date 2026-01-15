# Updates information for tiles after merging two tiles

Uses the information from the shared edge to update the merged PCs,
area, number of points, and perimeter. Mutates the `e_merge_from`th
value within the `V_pcs`, `V_npts`, `V_perimeter`, `V_area` input data
structures.

## Usage

``` r
update_V_cpp(
  V_pcs,
  V_npts,
  V_perimeter,
  V_area,
  e_merge_from,
  e_merge_to,
  e_merge_edge_length,
  e_merge_area,
  e_merge_npts,
  e_merge_pcs,
  agg_mode
)
```
