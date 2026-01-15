# Calculates triangles' centroids, areas, and heights from vertices

Calculates triangles' centroids, areas, and heights from vertices

## Usage

``` r
init_tris_cpp(triplets, pts)
```

## Arguments

- triplets:

  A M-by-3 matrix with indices for the points that correspond each
  triangle's vertices, where M is the number of triangles.

- pts:

  A N-by-2 matrix with indices for the X,Y coordinates of each point.

## Value

A M-by-4 matrix containing the X,Y coordinates of each triangle's
centroid in the first two columns, and area and largest height of each
triangle in the last two columns.
