# Compress a gradient field using SVD

Expresses the `2` x `D` total derivative at each location as a pair of
`2`-dimensional vectors in the gradient and orthogonal directions.

## Usage

``` r
compress_field_cpp(field)
```

## Arguments

- field:

  A `2` x `D` x `N` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every point, edge, or triangle.

## Value

A `N` x `6` matrix with the following attributes for each location:

- dx grad,dy grad:

  x,y directions of unit vector in the direction of greatest change
  (first singular vector).

- dx ortho,dy ortho:

  x,y directions of unit vector orthogonal to the direction of greatest
  change (second singular vector).

- \|grad\|,\|ortho\|:

  Magnitude of directional derivative in the gradient and orthogonal
  directions (singular values).
