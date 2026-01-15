# Compute spatial gradient field for input to DMT

First, spatial gradients for each embedding dimension are computed for
every points (cell) by looking at each cell's neighbors. These point
gradients can be smoothed using bilateral/anisotropic filtering. Then
gradients for triangles are computed as the average of their vertices (3
for proper triangles, 2 for degenerate exterior triangles at the
boundaries). Finally, primal and dual edge gradients are computed as the
sum of the gradients for the two points and two triangles, respectively,
that are associated with each edge.

## Usage

``` r
compute_gradients(
  dmt,
  smooth_distance = "none",
  smooth_similarity = "none",
  smooth_iter = 1,
  on_edges = FALSE
)
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

- smooth_distance:

  One of `c('none', 'euclidean', 'projected', 'constant')`. If either
  `smooth_distance` or `smooth_similarity` is `'none'` (the default),
  then no smoothing of the gradient field is conducted.

- smooth_similarity:

  One of `c('none', 'euclidean', 'projected', 'constant')`. If either
  `smooth_distance` or `smooth_similarity` is `'none'` (the default),
  then no smoothing of the gradient field is conducted.

- smooth_iter:

  Number of rounds of gradient smoothing.

## Value

A gradient field with the following attributes:

- pts:

  A `2` x `D` x `N` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every point in space.

- tris:

  A `2` x `D` x `F` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every triangle in the mesh. Average of the vertices (3 for full
  triangles, 2 for degenerate triangles).

- edges_pts:

  A `2` x `D` x `E` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every primal edge (point-to-point) in the mesh. Average of the two
  endpoints.

- edges_tris:

  A `2` x `D` x `E` array in column-major ordering containing the
  spatial gradient in expression for each of `D` latent variables at
  every dual edge (triangle-to-triangle) in the mesh. Average of the two
  adjacent triangles.
