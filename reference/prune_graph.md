# Prunes mesh by eliminating long edges and small connected components

Pruning works in three steps:

1.  Any triangle with an edge longer than a threshold length is removed.
    Afterwards, any edge that no longer belongs to a triangle is
    removed. Then any point that no longer belongs to an edge is
    removed.

2.  Connected components of triangles with a shared edge are computed.
    If everything is connected, then nothing is pruned.

3.  Otherwise, components (and corresponding triangles) that contain
    less than `mincells` points are removed. Afterwards, any edge that
    no longer belongs to a triangle is removed. Then any point that no
    longer belongs to an edge is removed.

## Usage

``` r
prune_graph(data, thresh_quantile = 0.95, mincells = 10, thresh = NA)
```

## Arguments

- data:

  A list containing the mesh data structures:

  - `pts` is a V-by-M data table with columns `X` and `Y` containing the
    coordinates of cells and additional metadata.

  - `tris` is a F-by-4 data table containing the X,Y coordinates of each
    triangle's centroid in the first two columns, and area and largest
    height of each triangle in the last two columns.

  - `edges` is a E-by-14 data table with columns `from_pt`, `to_pt`,
    `from_tri`, `to_tri`, `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`,
    `x1_tri`, `y0_tri`, `y1_tri`, `length_pt`, `length_tri`. If only one
    triangle uses an edge, then the `from_tri`, `x0_tri`, and `y0_tri`
    fields will contain NaN values.

  - `tri_to_pt` is a F-by-V sparse matrix with value 1 at (i,j) if
    triangle i uses point j as a vertex.

- thresh_quantile:

  Floating point value between 0 and 1, inclusive. Quantile of edge
  length above which edges are pruned. Defaults to 0.95.

- mincells:

  Minimum number of cells required for a connected component of
  triangles to be kept. Defaults to 10.

- thresh:

  Edge length above which edges are pruned. If equal to NA, then this
  value is ignored and `thresh_quantile` is used to compute the
  threshold. Otherwise, if `thresh` is set, then `thresh_quantile` is
  ignored. Defaults to NA.

## Value

A list containing the mesh data structures with (possibly) fewer points,
edges, and triangles. Indices have been updated since some objects might
have been removed.
