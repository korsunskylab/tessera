# Log-normalization for counts data

Normalizes all cells to have the same total counts and computes log(1+x)
transform.

## Usage

``` r
normalize_data(A, scaling_factor = NULL)
```

## Arguments

- A:

  A genes x cells counts matrix. Must be convertible to dgCMatrix.

- scaling_factor:

  Number of total counts to normalize all cells to. If NULL, then the
  median counts across all cells is used. Defaults to NULL.

## Value

Normalized genes x cells matrix in dgCMatrix sparse matrix format.
