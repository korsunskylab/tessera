# Computes the procrustes inner product between two matrices

Solves the orthogonal procrustes problem R = argmin \|\|R*A - B\|\|\_F
where R is an orthonormal matrix. Then computes the frobenius inner
product \<R*A, B\> = \<R, B\*A^T\>, which is maximized by R.

## Usage

``` r
procrustes_inner(A, B)
```

## Arguments

- A, B:

  Input matrices with same dimensions.

## Value

The procrustes inner product \<R\*A, B\> between A and B.
