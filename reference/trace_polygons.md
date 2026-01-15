# Contruct shapes that outline each tile

Contruct shapes that outline each tile

## Usage

``` r
trace_polygons(dmt, aggs)
```

## Arguments

- dmt:

  A DMT mesh data structure with `pts`, `edges`, and `tris` attributes.
  Assignment of points to tiles is provided in `dmt$pts$agg_id`.

- aggs:

  A tile data structure used to specify the total number of tiles.

## Value

A `sfc` list-column with the geometries for each tile.
