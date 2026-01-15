# Construct primal minimum spanning forest on points

Constructs a directed minimum spanning forest from point and edge scalar
values using a version of Prim's algorithm. Critical points (local
minimum; or, more precisely, an endpoint of an edge that is a local
minimum) are used as roots for each tree in the forest, and edges that
would bridge two trees with different critical point roots are marked as
possible saddle edges.

## Usage

``` r
do_primary_forest(dmt)
```

## Arguments

- dmt:

  A DMT data structure with the following attributes:

  - `pts$f`: Scalar field values defined for each point.

  - `edges$from_pt`: Index of first point joined by each edge.

  - `edges$to_pt`: Index of second point joined by each edge.

  - `edges$f_prim`: Scalar field values defined for each edge that
    connects two points.

## Value

A List with the following attributes (all indices are 1-indexed):

- `edges`: A data.table with `forest_size` rows, where each row is a
  directed edge in the minimum spanning forest. There are six columns:

  - `from,to`: Index of source and target points for each edge.

  - `x0,y0`: Coordinates of source point for each edge.

  - `x1,y1`: Coordinates of target point for each edge.

- `saddles`: A length `num_saddles` vector with edge indices for
  possible saddle edges.

- `labels`: A length `num_points` vector of labels for the connected
  components in the minimum spanning tree. Each connected component is
  labeled by the index of its critical point.

- `minima`: A length `num_critpts` vector of critical points (minima).

- `parent`: A length `num_points` vector containing the parent (source)
  point for each point in the directed spanning forest. Critical points
  have no parent, so the value should be ignored.

- `parent_edge`: A length `num_points` vector containing the directed
  edge that has each point as a target node. Critical points have no
  parent edge, so the value should be ignored.
