# Z-score a sparse matrix across each row

Z-score a sparse matrix across each row

## Usage

``` r
scaleRows_dgc(x, p, i, ncol, nrow, thresh)
```

## Arguments

- x, p, i:

  Internal data structures from a sparse matrix in dgCMatrix format.

- ncol, nrow:

  Dimensions of sparse matrix input.

- thresh:

  Z-scores above `thresh` and below `-thresh` are clipped to `thresh`
  and `-thresh`, respectively.

## Value

A dense matrix in column-major ordering with dimensions `nrow` x `ncol`.
