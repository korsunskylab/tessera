# Generate mesh data structure from coordinates, for input to DMT analysis

Creates a mesh from point coordinates using Delauney triangulation.
Stores the points, triangles, and edges for the mesh, as well as a
mapping from each triangle to its associated vertices. A gene-by-cell
matrix can also be included, as well as additional metadata for each
point.

## Usage

``` r
init_data(X, Y, counts, meta_data = NULL, meta_vars_include = c())
```

## Arguments

- X, Y:

  A pair of numeric vectors with the coordinates for each of V points.

- counts:

  A G x V gene-by-cell matrix of transcript counts.

- meta_data:

  A data frame with additional cell metadata to include.

- meta_vars_include:

  Names of columns in meta_data to include.

## Value

A list containing the mesh data structures:

- `pts` is a V-by-M data table with columns `X` and `Y` containing the
  coordinates of cells and additional metadata.

- `tris` is a F-by-4 data table containing the X,Y coordinates of each
  triangle's centroid in the first two columns, and area and largest
  height of each triangle in the last two columns.

- `edges` is a E-by-14 data table with columns `from_pt`, `to_pt`,
  `from_tri`, `to_tri`, `x0_pt`, `x1_pt`, `y0_pt`, `y1_pt`, `x0_tri`,
  `x1_tri`, `y0_tri`, `y1_tri`, `length_pt`, `length_tri`. If only one
  triangle uses an edge, then the `from_tri`, `x0_tri`, and `y0_tri`
  fields will contain NaN values.

- `tri_to_pt` is a F-by-V sparse matrix with value 1 at (i,j) if
  triangle i uses point j as a vertex.

- `counts` is a G x V gene-by-cell matrix of transcript counts.

Note that all indices stored in these data structures are 1-indexed.
