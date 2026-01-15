# Assign points to tiles after DMT

Tiles are defined as points that are in the same connected component in
the mesh after discounting separatrix edges.

## Usage

``` r
dmt_assign_tiles(dmt)
```

## Arguments

- dmt:

  A DMT data structure with `edges`, `pts`, and `e_sep` attributes.

## Value

A DMT data structure with the following additional attributes:

- `pts$agg_id`: Unique ID for the tile that each point belongs to.

- `edges$agg_from`: Unique ID for the tile that `edges$from_pt` belongs
  to.

- `edges$agg_to`: Unique ID for the tile that `edges$to_pt` belongs to.
