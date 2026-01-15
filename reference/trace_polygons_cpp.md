# Contruct shapes that outline each tile

Contruct shapes that outline each tile

## Usage

``` r
trace_polygons_cpp(edges, naggs, ntris, pts_dmt_component)
```

## Arguments

- edges:

  A matrix with `num_edges` rows with mesh edge information.

- naggs:

  Number of tiles.

- pts_dmt_component:

  A length `num_points` vector with the unique ID for the tile that each
  point belongs to.

- ntrix:

  Number of triangles in the mesh.

## Value

A list of length `naggs` matrices that contain coordinates of the
polygons that trace the outline of each tile.
