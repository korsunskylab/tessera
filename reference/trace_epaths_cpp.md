# Trace all paths from saddles to critical points in the spanning forest

Trace all paths from saddles to critical points in the spanning forest

## Usage

``` r
trace_epaths_cpp(saddles, vcrits, edges_from, edges_to, parent_edge, parent)
```

## Arguments

- saddles:

  A length `num_saddles` vector with edge indices for saddle edges
  (0-indexed).

- vcrits:

  A length `num_critpts` vector of critical points.

- edges_from:

  A length `num_edges` vector with indices for the first end-point of
  each edge in the mesh (0-indexed).

- edges_to:

  A length `num_edges` vector with indices for the second end-point of
  each edge in the mesh (0-indexed).

- parent_edge:

  A length `num_points` vector containing the directed edge that has
  each point as a target node in the directed spanning forest. Critical
  points have no parent edge, so the value is ignored.

- parent:

  A length `num_points` vector containing the parent (source) point for
  each point in the directed spanning forest. Critical points have no
  parent, so the value is ignored.

## Value

A list of length `2*num_saddles` containing the two paths from each
saddle edge to the two critical points that it joins. Each path is a
numeric vector of edge indices (1-indexed).
