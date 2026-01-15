# Update shapes and counts matrix for tiles after merging

Update shapes and counts matrix for tiles after merging

## Usage

``` r
update_agg_shapes(dmt, aggs)
```

## Arguments

- dmt:

  DMT data structure

- agg:

  Tile data structure after merging

## Value

`aggs` tile data structure with updated values for:

- meta_data\$shape:

  Retraces the polygons for tiles using the new point assignments.

- meta_data\$X,meta_data\$Y:

  Recomputes the centroid of each tile from the polygon shape.

- edges:

  Ensures that `aggs$edges` always has `from` \< `to`.

- counts:

  A `num_genes` x `num_tiles` gene-by-tile matrix of aggregated
  transcript counts.
