# Package index

## All functions

- [`AddAggsAdjacencyMatrix()`](https://korsunskylab.github.io/tessera/reference/AddAggsAdjacencyMatrix.md)
  : Construct tile adjacency matrix from consolidated GetTiles output.
- [`ConsolidateResults()`](https://korsunskylab.github.io/tessera/reference/ConsolidateResults.md)
  : Consolidate Tessera results from multiple samples (groups) after
  constructing Tessera tiles separately on cells from each group.
- [`GetTiles()`](https://korsunskylab.github.io/tessera/reference/GetTiles.md)
  : Generic function that runs the Tessera algorithm on single-cell
  spatial data
- [`GetTiles(`*`<Seurat>`*`)`](https://korsunskylab.github.io/tessera/reference/GetTiles.Seurat.md)
  : Applies Tessera on a Seurat object
- [`GetTiles(`*`<default>`*`)`](https://korsunskylab.github.io/tessera/reference/GetTiles.default.md)
  : Run full DMT segmentation pipeline to make aggregated tiles from
  cells
- [`RunUMAPCustom()`](https://korsunskylab.github.io/tessera/reference/RunUMAPCustom.md)
  : Run UMAP and save fgraph and embeddings in Seurat object
- [`add_exterior_triangles()`](https://korsunskylab.github.io/tessera/reference/add_exterior_triangles.md)
  : For every edge on the boundary, adds a second degenerate triangle
- [`assign_unique_rowid_cpp()`](https://korsunskylab.github.io/tessera/reference/assign_unique_rowid_cpp.md)
  : Assigns a unique ID to each point with distinct X,Y coordinates
- [`compress_field_cpp()`](https://korsunskylab.github.io/tessera/reference/compress_field_cpp.md)
  : Compress a gradient field using SVD
- [`compress_gradients_svd()`](https://korsunskylab.github.io/tessera/reference/compress_gradients_svd.md)
  : Compress a gradient field using SVD
- [`compute_gradients()`](https://korsunskylab.github.io/tessera/reference/compute_gradients.md)
  : Compute spatial gradient field for input to DMT
- [`compute_gradients_edges()`](https://korsunskylab.github.io/tessera/reference/compute_gradients_edges.md)
  : Compute spatial gradient field for input to DMT
- [`dmt_assign_tiles()`](https://korsunskylab.github.io/tessera/reference/dmt_assign_tiles.md)
  : Assign points to tiles after DMT
- [`dmt_get_separatrices()`](https://korsunskylab.github.io/tessera/reference/dmt_get_separatrices.md)
  : Get separatrices that separate points into components with strong
  boundaries
- [`dmt_init_tiles()`](https://korsunskylab.github.io/tessera/reference/dmt_init_tiles.md)
  : Initialize tiles with shapes and other properties
- [`dmt_set_f()`](https://korsunskylab.github.io/tessera/reference/dmt_set_f.md)
  : Set DMT scalar field values as the Frobenius norm of the total
  derivative
- [`do_dmt_forest_cpp()`](https://korsunskylab.github.io/tessera/reference/do_dmt_forest_cpp.md)
  : Construct maximum spanning forest
- [`do_dual_forest()`](https://korsunskylab.github.io/tessera/reference/do_dual_forest.md)
  : Construct dual maximum spanning forest on triangles
- [`do_pca()`](https://korsunskylab.github.io/tessera/reference/do_pca.md)
  : Compute PCA embeddings from a raw counts matrix
- [`do_primary_forest()`](https://korsunskylab.github.io/tessera/reference/do_primary_forest.md)
  : Construct primal minimum spanning forest on points
- [`estimate_field()`](https://korsunskylab.github.io/tessera/reference/estimate_field.md)
  : Compute a spatial gradient field at each point (cell)
- [`estimate_field_cpp()`](https://korsunskylab.github.io/tessera/reference/estimate_field_cpp.md)
  : Compute a spatial gradient field at each point (cell)
- [`estimate_field_edges_cpp()`](https://korsunskylab.github.io/tessera/reference/estimate_field_edges_cpp.md)
  : Compute a spatial gradient field along each edge
- [`findDuplicates()`](https://korsunskylab.github.io/tessera/reference/findDuplicates.md)
  : Find duplicates within a vector
- [`get_e_sep()`](https://korsunskylab.github.io/tessera/reference/get_e_sep.md)
  : Get the collection of edges that lie along separatrices
- [`init_data()`](https://korsunskylab.github.io/tessera/reference/init_data.md)
  : Generate mesh data structure from coordinates, for input to DMT
  analysis
- [`init_edges_cpp()`](https://korsunskylab.github.io/tessera/reference/init_edges_cpp.md)
  : Calculates triangles' centroids, areas, and heights from vertices
- [`init_scores()`](https://korsunskylab.github.io/tessera/reference/init_scores.md)
  : Initialize tile scores for aggregation
- [`init_tris_cpp()`](https://korsunskylab.github.io/tessera/reference/init_tris_cpp.md)
  : Calculates triangles' centroids, areas, and heights from vertices
- [`inv_svd_field_cpp()`](https://korsunskylab.github.io/tessera/reference/inv_svd_field_cpp.md)
  : Inverse SVD transform to reconstruct spatial gradient field
- [`mergeListsToArmaUVec()`](https://korsunskylab.github.io/tessera/reference/mergeListsToArmaUVec.md)
  : Copies elements from two lists an Armadillo uvec
- [`merge_aggs()`](https://korsunskylab.github.io/tessera/reference/merge_aggs.md)
  : Merges tiles using single-linkage agglomerative clustering
- [`merge_aggs_cpp()`](https://korsunskylab.github.io/tessera/reference/merge_aggs_cpp.md)
  : Merges tiles using single-linkage agglomerative clustering
- [`normalize_data()`](https://korsunskylab.github.io/tessera/reference/normalize_data.md)
  : Log-normalization for counts data
- [`procrustes_inner()`](https://korsunskylab.github.io/tessera/reference/procrustes_inner.md)
  : Computes the procrustes inner product between two matrices
- [`procrustes_mat()`](https://korsunskylab.github.io/tessera/reference/procrustes_mat.md)
  : Solves the orthogonal procrustes problem
- [`prune_graph()`](https://korsunskylab.github.io/tessera/reference/prune_graph.md)
  : Prunes mesh by eliminating long edges and small connected components
- [`scaleRows_dgc()`](https://korsunskylab.github.io/tessera/reference/scaleRows_dgc.md)
  : Z-score a sparse matrix across each row
- [`scale_data()`](https://korsunskylab.github.io/tessera/reference/scale_data.md)
  : Z-score a sparse matrix across each row or column
- [`smooth_embedding()`](https://korsunskylab.github.io/tessera/reference/smooth_embedding.md)
  : Smooth embeddings along edges
- [`smooth_field()`](https://korsunskylab.github.io/tessera/reference/smooth_field.md)
  : Bilateral / anisotropic filtering of gradient field
- [`smooth_field_cpp()`](https://korsunskylab.github.io/tessera/reference/smooth_field_cpp.md)
  : Bilateral / anisotropic filtering of gradient field
- [`smooth_field_edges()`](https://korsunskylab.github.io/tessera/reference/smooth_field_edges.md)
  : Bilateral / anisotropic filtering of gradient field
- [`smooth_field_edges_cpp()`](https://korsunskylab.github.io/tessera/reference/smooth_field_edges_cpp.md)
  : Bilateral / anisotropic filtering of gradient field
- [`svd_field_cpp()`](https://korsunskylab.github.io/tessera/reference/svd_field_cpp.md)
  : Compute the SVD of a spatial gradient field at each point (or
  edge/triangle)
- [`tessera-package`](https://korsunskylab.github.io/tessera/reference/tessera.md)
  [`tessera`](https://korsunskylab.github.io/tessera/reference/tessera.md)
  : Accurate tiling of spatial single-cell data
- [`tessera_warmup`](https://korsunskylab.github.io/tessera/reference/tessera_warmup.md)
  : Sample data for Tessera vignettes
- [`trace_back_cpp()`](https://korsunskylab.github.io/tessera/reference/trace_back_cpp.md)
  : Trace back from a point to its root in the spanning forest
- [`trace_epaths_cpp()`](https://korsunskylab.github.io/tessera/reference/trace_epaths_cpp.md)
  : Trace all paths from saddles to critical points in the spanning
  forest
- [`trace_paths()`](https://korsunskylab.github.io/tessera/reference/trace_paths.md)
  : Trace all paths from saddles to dual critical points in the spanning
  forest.
- [`trace_polygons()`](https://korsunskylab.github.io/tessera/reference/trace_polygons.md)
  : Contruct shapes that outline each tile
- [`trace_polygons_cpp()`](https://korsunskylab.github.io/tessera/reference/trace_polygons_cpp.md)
  : Contruct shapes that outline each tile
- [`update_E_cpp()`](https://korsunskylab.github.io/tessera/reference/update_E_cpp.md)
  : Updates information for boundaries after merging two tiles
- [`update_V_cpp()`](https://korsunskylab.github.io/tessera/reference/update_V_cpp.md)
  : Updates information for tiles after merging two tiles
- [`update_agg_shapes()`](https://korsunskylab.github.io/tessera/reference/update_agg_shapes.md)
  : Update shapes and counts matrix for tiles after merging
- [`update_dmt_aggid()`](https://korsunskylab.github.io/tessera/reference/update_dmt_aggid.md)
  : Update tile IDs in the DMT data structure after merging
