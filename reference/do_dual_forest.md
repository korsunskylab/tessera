# Construct dual maximum spanning forest on triangles

Constructs a directed maximum spanning forest from triangle and edge
scalar values using a version of Prim's algorithm. Critical triangles
(local maximum; or, more precisely, a triangle that contains an edge
that is a local maximum) are used as roots for each tree in the forest,
and edges that would bridge two trees with different critical triangle
roots are marked as possible saddle edges.

## Usage

``` r
do_dual_forest(dmt)
```

## Arguments

- dmt:

  A DMT data structure with the following attributes:

  - `tris$f`: Scalar field values defined for each triangle.

  - `edges$from_tri`: Index of first triangle joined by each edge.

  - `edges$to_tri`: Index of second triangle joined by each edge.

  - `edges$f_dual`: Scalar field values defined for each edge that
    connects two triangles.

## Value

A List with the following attributes (all indices are 1-indexed):

- `edges`: A data.table with `forest_size` rows, where each row is a
  directed edge in the maximum spanning forest. There are six columns:

  - `from,to`: Index of source and target triangles for each edge.

  - `x0,y0`: Coordinates of source triangle for each edge.

  - `x1,y1`: Coordinates of target triangle for each edge.

- `saddles`: A length `num_saddles` vector with edge indices for
  possible saddle edges.

- `labels`: A length `num_triangles` vector of labels for the connected
  components in the maximum spanning tree. Each connected component is
  labeled by the index of its critical triangle.

- `maxima`: A length `num_critpts` vector of critical triangles
  (maxima).

- `parent`: A length `num_triangles` vector containing the parent
  (source) triangle for each triangle in the directed spanning forest.
  Critical triangles have no parent, so the value should be ignored.

- `parent_edge`: A length `num_triangles` vector containing the directed
  edge that has each triangle as a target node. Critical triangles have
  no parent edge, so the value should be ignored.
