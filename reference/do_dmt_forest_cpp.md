# Construct maximum spanning forest

Constructs a directed maximum spanning forest from point and edge scalar
values using a version of Prim's algorithm. Critical points (local
maxima; or, more precisely, an endpoint of an edge that is a local
maxima) are used as roots for each tree in the forest, and edges that
would bridge two trees with different critical point roots are marked as
possible saddle edges.

## Usage

``` r
do_dmt_forest_cpp(f, edges_from, edges_to, edges_f)
```

## Arguments

- f:

  Vector of `num_points` scalar values defined at each point.

- edges_from:

  Vector of `num_edges` indices for first end-point of each edge
  (0-indexed).

- edges_to:

  Vector of `num_edges` indices for second end-point of each edge
  (0-indexed).

- edges_f:

  Vector of `num_edges` scalar values defined along each edge.

## Value

A List with the following attributes (all indices are 1-indexed):

- `edges`: A `forest_size` x `2` matrix where each row is a directed
  edge in the maximum spanning forest. The first column has the source
  point for each edge and the second column has the target point.

- `saddles`: A length `num_saddles` vector with edge indices for
  possible saddle edges.

- `labels`: A length `num_points` vector of labels for the connected
  components in the maximum spanning tree. Each connected component is
  labeled by the index of its critical point.

- `critpts`: A length `num_critpts` vector of critical points (maxima).

- `parent`: A length `num_points` vector containing the parent (source)
  point for each point in the directed spanning forest. Critical points
  have no parent, so the value should be ignored.

- `parent_edge`: A length `num_points` vector containing the directed
  edge that has each point as a target node. Critical points have no
  parent edge, so the value should be ignored.
