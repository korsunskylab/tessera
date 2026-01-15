# Inverse SVD transform to reconstruct spatial gradient field

Inverse SVD transform to reconstruct spatial gradient field

## Usage

``` r
inv_svd_field_cpp(svd_list)
```

## Arguments

- svd_list:

  A list with elements `u` (`2` x `2` x `N`), `s` (`2` x `N`), and `v`
  (`D` x `2` x `N`) containing the left singular vectors, singular
  values, and right singular vectors for each point (or edge/triangle).

## Value

A `2` x `D` x `N` array in column-major ordering containing the spatial
gradient in expression for each of `D` embedding dimensions at every
point (or edge/triangle).
