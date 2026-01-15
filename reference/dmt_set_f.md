# Set DMT scalar field values as the Frobenius norm of the total derivative

Set DMT scalar field values as the Frobenius norm of the total
derivative

## Usage

``` r
dmt_set_f(dmt, field)
```

## Arguments

- dmt:

  A list containing the mesh data structures:

  - `pts` is a `N` x `2+M` data table with columns `X` and `Y`
    containing the coordinates of cells and additional metadata.

  - `tris` is a `F` x `4` data table containing the X,Y coordinates of
    each triangle's centroid in the first two columns, and area and
    largest height of each triangle in the last two columns.

  - `edges` is a `E` x `14` data table with columns `from_pt`, `to_pt`,
    `from_tri`, `to_tri`, `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`,
    `x1_tri`, `y0_tri`, `y1_tri`, `length_pt`, `length_tri`. If only one
    triangle uses an edge, then the `from_tri`, `x0_tri`, and `y0_tri`
    fields will contain NaN values.

  - `tri_to_pt` is a `F` x `N` sparse matrix with value 1 at (i,j) if
    triangle i uses point j as a vertex.

  - `udv_cells` contains cell embeddings stored in `embeddings` (a `N` x
    `D` matrix with `D`-dimensional embeddings for each cell) and
    `loadings` (a `G` x `D` matrix with gene loadings for each latent
    variable).

- field:

  A gradient field with the compressed representations `pts_svd`,
  `tris_svd`, `edges_pts_svd`, and `edges_tris_svd`. Each of these is a
  `N` x `6` matrix with the following columns for each location:

  - `dx_grad,dy_grad`: x,y directions of unit vector in the direction of
    greatest change (first singular vector).

  - `dx_ortho,dy_ortho`: x,y directions of unit vector orthogonal to the
    direction of greatest change (second singular vector).

  - `len_grad,len_ortho`: Magnitude of directional derivative in the
    gradient and orthogonal directions (singular values).

## Value

`dmt` with the following additional attributes: `dmt$pts$f`,
`dmt$tris$f`, `dmt$edges$f_prim`, and `dmt$edges$f_dual`
