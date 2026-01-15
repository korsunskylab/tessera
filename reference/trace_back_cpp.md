# Trace back from a point to its root in the spanning forest

Trace back from a point to its root in the spanning forest

## Usage

``` r
trace_back_cpp(v0, vcrit, parent_edge, parent)
```

## Arguments

- v0:

  Index of starting point (0-indexed).

- vcrit:

  Index of critical point associated with the tree that `v0` belongs to
  (0-indexed).

- parent_edge:

  A length `num_points` vector containing the directed edge that has
  each point as a target node. Critical points have no parent edge, so
  the value is ignored.

- parent:

  A length `num_points` vector containing the parent (source) point for
  each point in the directed spanning forest. Critical points have no
  parent, so the value is ignored.

## Value

Vector of edge indices along the path from `v0` to its root `vcrit` in
the tree (1-indexed).
