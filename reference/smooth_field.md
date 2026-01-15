# Bilateral / anisotropic filtering of gradient field

Gradient fields are smoothed using bilateral filtering, in which the
smoothed gradient of each point is computed as the weighted average of
the neighbors' gradients, considering both distance in space and also
similarity in gradients.

## Usage

``` r
smooth_field(
  coords,
  field,
  adj,
  include_self = TRUE,
  distance = "euclidean",
  similarity = "euclidean"
)
```

## Arguments

- coords:

  A `N` x `2` matrix of cell coordinates.

- field:

  A `2` x `D` x `N` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every point in space.

- adj:

  A `N` x `N` sparse adjacency matrix in dgCMatrix format.

- include_self:

  A boolean whether or not to include the each point's gradient in its
  own smoothed value. Defaults to TRUE.

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

## Details

The weight of each neighbor is computed from the product of two scores:

- `distance` score: Generally, closer neighbors have greater weight.

  - if `'euclidean'`: Gaussian transformation of the Euclidean distance
    of a cell from its neighbor, so that more distant neighbors have
    less weight.

  - if `'projected'`: An anisotropic filter that accounts for expected
    change in expression along the direction of the neighbor. The
    expected change in expression is calculated from the gradient field
    as the total derivative in the direction of the neighbor. This
    change in expression is then Gaussian transformed so that neighbors
    that are more distant along the direction of greatest change have
    less weight.

  - if `'constant'`: All neighbors have equal `distance` weights

- `similarity` score: Generally, neighbors with more similar gradients
  have greater weight

  - if `'euclidean'`: Gaussian transformation of the Euclidean distance
    between a cell's gradient field and its neighbor's gradient field.

  - if `'projected'`: Gaussian transformation of the cosine distance
    between a cell's gradient field and its neighbor's gradient field.

  - if `'constant'`: All neighbors have equal `similarity` weights
