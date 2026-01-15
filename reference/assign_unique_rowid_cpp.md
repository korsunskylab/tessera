# Assigns a unique ID to each point with distinct X,Y coordinates

Assigns a unique ID to each point with distinct X,Y coordinates

## Usage

``` r
assign_unique_rowid_cpp(X, Y)
```

## Arguments

- X, Y:

  A pair of numeric vectors with the coordinates for each point.

## Value

A vector with the same length as `X` and `Y` containing IDs that range
from 0 to N where N is the number of unique points.
