# Get the collection of edges that lie along separatrices

Get the collection of edges that lie along separatrices

## Usage

``` r
get_e_sep(epaths, saddles, nedges)
```

## Arguments

- epaths:

  A list of length `2*num_saddles` containing the two paths from each
  saddle edge to the two critical points that it joins. Each path is a
  numeric vector of edge indices (1-indexed).

- saddles:

  A length `num_saddles` vector with edge indices for saddle edges
  (1-indexed).

- nedges:

  Total number of edges in the mesh.

## Value

A length `num_sep_edges` vector of edges (0-indexed) that are saddle
edges in `saddles` or lie along the paths in `epaths`. These edges make
up the separatrices that separate points into different components.
