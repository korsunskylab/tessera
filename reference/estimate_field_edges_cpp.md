# Compute a spatial gradient field along each edge

Distance between neighboring cells is normalized to unit distance so
that only the orientation of each edge matters. The gradient is then the
difference in expression of each embedding dimension between the two
endpoints of the edge.

## Usage

``` r
estimate_field_edges_cpp(coords, embeddings, from_pt, to_pt)
```

## Arguments

- coords:

  A `N` x `2` matrix of cell coordinates.

- embeddings:

  A `N` x `D` matrix of cell embeddings.

- from_pt, to_pt:

  A pair of `E`-length vectors indicating the start and end points of
  each edge (0-indexed).

## Value

A `2` x `D` x `E` array in column-major ordering containing the spatial
gradient in expression for each of `D` embedding dimensions at every
edge.
