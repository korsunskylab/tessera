# Update tile IDs in the DMT data structure after merging

Update tile IDs in the DMT data structure after merging

## Usage

``` r
update_dmt_aggid(dmt, aggs)
```

## Arguments

- dmt:

  DMT data structure

- agg:

  Tile data structure after merging, with `aggmap` attribute mapping
  original tile IDs to new IDs for the merged tiles.

## Value

DMT data structure with updated values for these attributes:

- `pts$agg_id`: Unique ID for the tile that each point belongs to.

- `edges$agg_from`: Unique ID for the tile that `edges$from_pt` belongs
  to.

- `edges$agg_to`: Unique ID for the tile that `edges$to_pt` belongs to.
