# Smooth embeddings along edges

Smooth embeddings along edges

## Usage

``` r
smooth_embedding(dmt, smooth_emb = 0)
```

## Arguments

- dmt:

  A DMT object.

- smooth_emb:

  Number of smoothing iterations to perform. If a vector, then
  embeddings after each specified iteration are concatenated. If `0` is
  included, then the original embeddings are also included.

## Value

The input DMT object with smoothed embeddings stored in
`dmt$udv_cells$embeddings`.
