# Get separatrices that separate points into components with strong boundaries

Finds the the collection of edges that lie along separatrices. These are
found by tracing paths from saddles to critical triangles (maxima) in
the dual spanning forest. Because the dual spanning forest is a maximum
spanning forest on triangles in the mesh, these paths will follow
ridges, separating points into components with the strongest boundaries.

## Usage

``` r
dmt_get_separatrices(dmt)
```

## Arguments

- dmt:

  A DMT data structure with the following attributes:

  - `edges`: Data structure of edges in the mesh.

  - `prim`: Data structure of primal minimum spanning forest.

  - `dual`: Data structure of dual maximum spanning forest.

## Value

A length `num_sep_edges` vector of edges (1-indexed) that make up the
separatrices that separate points into different components.

## Details

Saddle edges are found by interesting the saddle edges found in the
primal forest with saddles edges found in the dual forest. All boundary
edges in the dual forest are also included as possible saddle edges.
