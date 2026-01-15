# Trace all paths from saddles to dual critical points in the spanning forest.

Because the dual spanning forest is a maximum spanning forest on
triangles in the mesh, the paths that are traced between saddles and
critical triangles (maxima) will follow ridges, separating points into
components with the strongest boundaries.

## Usage

``` r
trace_paths(dmt, dual, saddles)
```

## Arguments

- dmt:

  A DMT data structure with the following attributes:

  - `edges$from_tri`: Index of first triangle joined by each edge
    (1-indexed).

  - `edges$to_tri`: Index of second triangle joined by each edge
    (1-indexed).

- dual:

  A dual maximum spanning forest on triangles:

  - `labels`: A length `num_triangles` vector of labels for the
    connected components in the maximum spanning tree. Each connected
    component is labeled by the index of its critical triangle.

  - `parent`: A length `num_triangles` vector containing the parent
    (source) triangle for each triangle in the directed spanning forest.
    Critical triangles have no parent, so the value should be ignored.

  - `parent_edge`: A length `num_triangles` vector containing the
    directed edge that has each triangle as a target node. Critical
    triangles have no parent edge, so the value should be ignored.

- saddles:

  A length `num_saddles` vector with edge indices for saddle edges
  (1-indexed).

## Value

A list of length `2*num_saddles` containing the two paths from each
saddle edge to the two critical points that it joins. Each path is a
numeric vector of edge indices (1-indexed).
