# Z-score a sparse matrix across each row or column

Z-score a sparse matrix across each row or column

## Usage

``` r
scale_data(A, margin = 1, thresh = 10)
```

## Arguments

- A:

  A log-transformed genes x cells counts matrix. Must be convertible to
  dgCMatrix.

- margin:

  Subscript over which to compute Z-score. If `1`, then each row is
  Z-scored; otherwise, each column is Z-scored.

## Value

Z-scored genes x cells matrix in dgCMatrix sparse matrix format.
