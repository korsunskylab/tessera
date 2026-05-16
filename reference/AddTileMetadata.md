# Add tile-level metadata to a matching cell-level Seurat object

Add tile-level metadata to a matching cell-level Seurat object

## Usage

``` r
AddTileMetadata(obj, tile_obj, ..., by = c(tile_id = "id"))
```

## Arguments

- obj:

  A Seurat object containing cell-level data with a `tile_id` metadata
  column.

- tile_obj:

  A Seurat object containing tile-level data.

- ...:

  Tidy expressions specifying which metadata columns from `tile_obj` to
  add to `obj`. For example, `tile_size = npts`

- by:

  A named character vector specifying how to join `obj` and `tile_obj`.
  Defaults to `c('tile_id' = 'id')`, where `tile_id` is the metadata
  column in `obj` and `id` is the metadata column in `tile_obj`.

## Value

The input Seurat object with additional tile-level metadata added to
`@meta.data`.
