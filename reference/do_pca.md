# Compute PCA embeddings from a raw counts matrix

First, log-normalizes and Z-scores the counts matrix and then performs
PCA using SVD.

## Usage

``` r
do_pca(counts, npcs)
```

## Arguments

- counts:

  A `n_genes` x `n_cells` counts matrix. Must be convertible to
  dgCMatrix.

- npcs:

  Number of PCs to compute.

## Value

A list with two features:

- `loadings`: A `n_genes` x `npcs` matrix of gene loadings for each PC.
  Each column is a unit vector.

- `embeddings`: A `n_cells` x `npcs` matrix of cell embeddings across
  all PCs. Each column *j* has magnitude equal to the *j*th singular
  value. That is, PCs with larger contribution to the total variance
  will have embeddings of proportionally larger magnitude.
