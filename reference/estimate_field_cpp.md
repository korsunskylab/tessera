# Compute a spatial gradient field at each point (cell)

Distance between neighboring cells is normalized to unit distance so
that only the direction from each cell to its neighbors matters. The
gradient is then the average gradient in expression of each embedding
dimension between the index cell and its neighbors.

## Usage

``` r
estimate_field_cpp(coords, embeddings, adj_i, adj_p)
```

## Arguments

- coords:

  A `N` x `2` matrix of cell coordinates.

- embeddings:

  A `N` x `D` matrix of cell embeddings.

- adj_i, adj_p:

  A `N` x `N` sparse adjacency matrix in dgCMatrix format.

## Value

A `2` x `D` x `N` array in column-major ordering containing the spatial
gradient in expression for each of `D` embedding dimensions at every
point in space.
