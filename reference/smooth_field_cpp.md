# Bilateral / anisotropic filtering of gradient field

Gradient fields are smoothed using bilateral filtering, in which the
smoothed gradient of each point is computed as the weighted average of
the neighbors' gradients, considering both distance in space and also
similarity in gradients.

## Usage

``` r
smooth_field_cpp(pvec, adj_i, adj_p, field, coords, distance, similarity)
```

## Arguments

- pvec, adj_i, adj_p:

  A `N` x `N` sparse adjacency matrix in dgCMatrix format:
  `pvec = diff(adj@p)`, `adj_i = adj@i`, and `adj_p = adj@p`

- field:

  A `2` x `D` x `N` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every point in space.

- coords:

  A `N` x `2` matrix of cell coordinates.

- distance:

  Method for computing distance score in weighted average. See
  description for details. Defaults to `'euclidean'`.

- similarity:

  Method for computing similarity score in weighted average. See
  description for details. Defaults to `'euclidean'`.

## Value

A `2` x `D` x `N` array in column-major ordering containing the smoothed
spatial gradient in expression for each of `D` latent variables at every
point in space.
