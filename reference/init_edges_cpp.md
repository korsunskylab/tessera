# Calculates triangles' centroids, areas, and heights from vertices

Calculates triangles' centroids, areas, and heights from vertices

## Usage

``` r
init_edges_cpp(triplets, pts, tris)
```

## Arguments

- triplets:

  A M-by-3 matrix with indices for the points that correspond each
  triangle's vertices, where M is the number of triangles.

- pts:

  A N-by-2 matrix with indices for the X,Y coordinates of each point.

- tris:

  A M-by-4 matrix containing the X,Y coordinates of each triangle's
  centroid in the first two columns, and area and largest height of each
  triangle in the last two columns.

## Value

A E-by-14 matrix with columns `from_pt`, `to_pt`, `from_tri`, `to_tri`,
`x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`, `x1_tri`, `y0_tri`,
`y1_tri`, `length_pt`, `length_tri`. If only one triangle uses an edge,
then the `from_tri`, `x0_tri`, and `y0_tri` fields will contain NaN
values.

Note that indices that reference the `pts` and `tris` tables are
1-indexed.
