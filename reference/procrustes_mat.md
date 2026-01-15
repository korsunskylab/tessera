# Solves the orthogonal procrustes problem

Solves the orthogonal procrustes problem R = argmin \|\|R\*A - B\|\|\_F
where R is an orthonormal matrix.

## Usage

``` r
procrustes_mat(A, B)
```

## Arguments

- A, B:

  Input matrices with same dimensions.

## Value

The orthonormal matrix R that best maps A to B.
