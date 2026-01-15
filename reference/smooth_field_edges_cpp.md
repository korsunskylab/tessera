# Bilateral / anisotropic filtering of gradient field

Gradient fields are smoothed using bilateral filtering, in which the
smoothed gradient of each edge is computed as the weighted average of
the neighboring edges' gradients, considering both distance in space and
also similarity in gradients.

## Usage

``` r
smooth_field_edges_cpp(
  from_pt,
  to_pt,
  field,
  edges_svd,
  coords,
  adj_idx,
  distance,
  similarity
)
```

## Arguments

- from_pt, to_pt:

  A pair of `E`-length vectors indicating the start and end points of
  each edge (0-indexed).

- field:

  A `2` x `D` x `E` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every edge.

- edges_svd:

  A list with elements `u` (`2` x `2` x `E`), `s` (`2` x `E`), and `v`
  (`D` x `2` x `E`) containing the left singular vectors, singular
  values, and right singular vectors for each edge.

- coords:

  A `N` x `2` matrix of cell coordinates.

- adj_idx:

  A `N` x `N` sparse adjacency matrix in dgCMatrix format, containing
  the mapping from cells to edges. In particular, `adj_idx@x` should
  store the 1-indexed edge indices corresponding to each cell-cell pair.
  The adjacency matrix is assumed to be symmetric and have zeros on the
  diagonal. Importantly, note that the stored edge indices are 1-indexed
  (*not* 0-indexed here) in order to avoid problems with R's sparse
  matrix representation.

- distance:

  Method for computing distance score in weighted average. See
  description for details. Defaults to `'euclidean'`.

- similarity:

  Method for computing similarity score in weighted average. See
  description for details. Defaults to `'euclidean'`.

## Value

A `2` x `D` x `E` array in column-major ordering containing the smoothed
spatial gradient in expression for each of `D` latent variables at every
edge.
