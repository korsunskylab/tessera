# For every edge on the boundary, adds a second degenerate triangle

Edges at the boundary will only border a single triangle. In order to
deal with boundary conditions properly, this step ensures that every
edge has two associated triangles (which use the edge) as well as two
associated points (endpoints). The triangles that are added are
degenerate triangles centered at the midpoint of each boundary edge.

## Usage

``` r
add_exterior_triangles(data)
```

## Arguments

- data:

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

## Value

A list containing the mesh data structures with (possibly) additional
triangles. The `tris` table is updated to include the added external
triangles. The `edge` table is updated so that boundary edges point to
the newly added triangles. Also `pts` is unchanged, `tri_to_pt` is
updated. Because new external triangles only have 2 associated points,
`tri_to_pt` has rows that sum to either 2 or 3.
