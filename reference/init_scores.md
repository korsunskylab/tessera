# Initialize tile scores for aggregation

Higher scores favor merging.

## Usage

``` r
init_scores(aggs, agg_mode, ...)
```

## Arguments

- aggs:

  A tile data structure.

- agg_mode:

  Method to use to calculate aggregation scores (1, 2, or 3).

- ...:

  Additional parameters for different modes:

  1.  No additional parameters.

  2.  Requires `alpha` and `max_npts`. For `alpha`, 0.2 = conservative
      merging, 2 = liberal merging.

  3.  Requires `alpha` and `min_npts`. For `alpha`, 0.2 = conservative
      merging, 2 = liberal merging.

## Value

`aggs` with scores for aggregation stored in attributes:

- edges:

  Additional attributes are calculated: \* `dscore`: Overall score for
  merging two tiles. Product of `w`, `score_size`, and `dC`. \* `w`:
  Gene expression similarity score. \* `score_size`: Penalizes tiles
  with many points. \* `perimeter_merge`: Perimeter of merged tile.

- d_mu,d_sig:

  Parameters used to calculate `w`.

- pcs_merged:

  Average PCs for merged tile.

## Details

Methods for different modes of aggregation:

1.  `dscore` has already been manually calculated and stored in
    `aggs$edges$dscore`. Will set the `score`, `score_size`, and
    `compactness` attributes of `aggs$meta_data` to 0.

2.  `dscore` is the product of three factors:

    - `w`: A 2-cluster GMM is used to determine the mean `mu` and
      standard deviation `sig` of the distance between tiles that have
      similar gene expression (Euclidean distance `d` in PC space). Then
      we define `d_mu = mu + sig` and `d_sig = alpha * sig`, and
      calculate `w` as `w = 0.5 - 1 / (1 + exp(-(d - d_mu) / d_sig))`.
      Ranges from -0.5 to 0.5. If adjacent tiles are very dissimilar
      (`d >> d_mu`), then `d` is large, and `w` is close to `-0.5`. If
      adjacent tiles are very similar (`d < d_mu`), then `d` is small,
      and `w` is positive.

    - `score_size`: `(1 - npts_from/max_npts) * (1 - npts_to/max_npts)`.
      Ranges from 0 to 1. If merging the two tiles would have a total
      number of points ≥`max_npts`, then `score_size` is `-Inf`, which
      prevents merging.

    - `dC`: `.5 * (C_merge - C_from - C_to + 1)`. Ranges from 0 to 1.

3.  `dscore` is the product of same three factors as mode 2, but
    `score_size` has already been precomputed. Additionally, `dscore` is
    set to -1 if both adjacent tiles have at least `min_npts` cells.
