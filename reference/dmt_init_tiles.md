# Initialize tiles with shapes and other properties

Initialize tiles with shapes and other properties

## Usage

``` r
dmt_init_tiles(dmt)
```

## Arguments

- dmt:

  A DMT data structure with `pts`, `edges`, and `udv_cells` attributes.
  Assignment of points to tiles is provided in `dmt$pts$agg_id`,
  `dmt$edges$agg_from`, and `dmt$edges$agg_to`.

## Value

A list data structure for tiles and their adjacencies. Includes the
following attributes:

- `meta_data`: A data table with attributes for tiles:

  - `X,Y`: Centroid of each tile.

  - `npts`: Number of points in each tile.

  - `shape`: A `sfc` list-column with the geometries for each tile.

  - `area,perimeter`: Area and perimeter of each tile.

- `edges`: A data table with attributes for edges between adjacent
  tiles:

  - `from,to`: Tile IDs for the two tiles bordering this edge.

  - `x0,y0,x1,y1`: Centroid coordinates for the two tiles bordering this
    edge.

  - `area,npts`: Sum of areas and numbers of points in the two tiles
    bordering this edge.

  - `edge_length`: Total length of the border between the `from` and
    `to` tiles.

- `pcs`: A `num_tiles` x `npcs` matrix with the average embedding value
  over all cells in each tile.
