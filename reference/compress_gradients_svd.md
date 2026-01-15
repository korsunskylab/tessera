# Compress a gradient field using SVD

Expresses the `2` x `D` total derivative at each location as a pair of
`2`-dimensional vectors in the gradient and orthogonal directions.

## Usage

``` r
compress_gradients_svd(field)
```

## Arguments

- field:

  A gradient field with the following attributes:

  - `pts`: A `2` x `D` x `N` array in column-major ordering containing
    the spatial gradient in expression for each of `D` latent variables
    at every point in space.

  - `tris`: A `2` x `D` x `F` array in column-major ordering containing
    the spatial gradient in expression for each of `D` latent variables
    at every triangle in the mesh. Average of the vertices (3 for full
    triangles, 2 for degenerate triangles).

  - `edges_pts`: A `2` x `D` x `E` array in column-major ordering
    containing the spatial gradient in expression for each of `D` latent
    variables at every primal edge (point-to-point) in the mesh. Sum of
    the two endpoints.

  - `edges_tris`: A `2` x `D` x `E` array in column-major ordering
    containing the spatial gradient in expression for each of `D` latent
    variables at every dual edge (triangle-to-triangle) in the mesh. Sum
    of the two adjacent triangles.

## Value

A gradient field with the same attributes as the input, as well as
compressed representations `pts_svd`, `tris_svd`, `edges_pts_svd`, and
`edges_tris_svd`. Each of these is a `N` x `6` matrix with the following
columns for each location:

- dx_grad,dy_grad:

  x,y directions of unit vector in the direction of greatest change
  (first singular vector).

- dx_ortho,dy_ortho:

  x,y directions of unit vector orthogonal to the direction of greatest
  change (second singular vector).

- len_grad,len_ortho:

  Magnitude of directional derivative in the gradient and orthogonal
  directions (singular values).
