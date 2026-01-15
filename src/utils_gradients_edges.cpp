#include <RcppArmadillo.h>
using namespace Rcpp;

//' Compute a spatial gradient field along each edge
//'
//' Distance between neighboring cells is normalized to unit
//' distance so that only the orientation of each edge matters.
//' The gradient is then the difference in expression of each
//' embedding dimension between the two endpoints of the edge.
//'
//' @param coords A `N` x `2` matrix of cell coordinates.
//' @param embeddings A `N` x `D` matrix of cell embeddings.
//' @param from_pt,to_pt A pair of `E`-length vectors indicating
//'   the start and end points of each edge (0-indexed).
//'
//' @returns A `2` x `D` x `E` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` embedding dimensions at every edge.
//'
// [[Rcpp::export]]
arma::cube estimate_field_edges_cpp(
    arma::mat & coords, 
    arma::mat & embeddings, 
    arma::uvec & from_pt, 
    arma::uvec & to_pt
) {
    unsigned N = from_pt.n_elem; 
    unsigned D = embeddings.n_cols; 
    arma::cube field = arma::zeros<arma::cube>(2, D, N); 

    // compute embeddings differences along edges
    arma::mat f = embeddings.rows(to_pt) - embeddings.rows(from_pt);

    // center and scale the XY coordinates
    arma::mat xy = coords.rows(to_pt) - coords.rows(from_pt);
    xy = arma::normalise(xy, 2, 1); // p=2, dim=1 (rows)

    for (unsigned i = 0; i < N; i++) {
        field.slice(i) = xy.row(i).t() * f.row(i);
    } 
        
    return field; 
}

//' Compute the SVD of a spatial gradient field at each point (or edge/triangle)
//'
//' @param field A `2` x `D` x `N` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` embedding dimensions at every point (or edge/triangle).
//'
//' @returns A list with elements `u` (`2` x `2` x `N`), `s` (`2` x `N`), and
//'   `v` (`D` x `2` x `N`) containing the left singular vectors, singular values,
//'   and right singular vectors for each point (or edge/triangle).
//'
// [[Rcpp::export]]
Rcpp::List svd_field_cpp(
    arma::cube & field    
) {
    unsigned D = field.n_cols; 
    unsigned N = field.n_slices; 
    arma::cube u = arma::zeros<arma::cube>(2, 2, N);
    arma::mat s = arma::zeros<arma::mat>(2, N);
    arma::cube v = arma::zeros<arma::cube>(D, 2, N);

    arma::mat u_i, v_i;
    arma::vec s_i;
    for (unsigned i = 0; i < N; i++) {
        arma::svd_econ(u_i, s_i, v_i, field.slice(i)); 
        u.slice(i) = u_i;
        s.col(i) = s_i;
        v.slice(i) = v_i;
    }
    return Rcpp::List::create(
        _["u"]  = u,
        _["s"]  = s,
        _["v"]  = v
    );
}

//' Inverse SVD transform to reconstruct spatial gradient field
//'
//' @param svd_list A list with elements `u` (`2` x `2` x `N`), `s` (`2` x `N`), and
//'   `v` (`D` x `2` x `N`) containing the left singular vectors, singular values,
//'   and right singular vectors for each point (or edge/triangle).
//'
//' @returns A `2` x `D` x `N` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` embedding dimensions at every point (or edge/triangle).
//'
// [[Rcpp::export]]
arma::cube inv_svd_field_cpp(
    Rcpp::List svd_list
) {
    arma::cube u = svd_list["u"]; // (2, 2, N)
    arma::mat s = svd_list["s"];  // (2, N)
    arma::cube v = svd_list["v"]; // (D, 2, N)
    
    unsigned D = v.n_rows; 
    unsigned N = v.n_slices; 
    arma::cube field = arma::zeros<arma::cube>(2, D, N);

    for (unsigned i = 0; i < N; i++) {
        field.slice(i) = u.slice(i) * diagmat(s.col(i)) * v.slice(i).t();
    }
    return field;
}

//' Solves the orthogonal procrustes problem
//'
//' Solves the orthogonal procrustes problem
//'     R = argmin ||R*A - B||_F
//' where R is an orthonormal matrix.
//'
//' @param A,B Input matrices with same dimensions.
//'
//' @returns The orthonormal matrix R that best maps A to B.
//'
// [[Rcpp::export]]
arma::mat procrustes_mat(
    arma::mat A,
    arma::mat B
) {
    arma::mat u, v;
    arma::vec s;
    arma::svd_econ(u, s, v, B * A.t());

    return u * v.t();
}

//' Computes the procrustes inner product between two matrices
//' 
//' Solves the orthogonal procrustes problem
//'     R = argmin ||R*A - B||_F
//' where R is an orthonormal matrix.
//' Then computes the frobenius inner product
//'     <R*A, B> = <R, B*A^T>,
//' which is maximized by R.
//'
//' @param A,B Input matrices with same dimensions.
//'
//' @returns The procrustes inner product <R*A, B> between A and B.
//'
// [[Rcpp::export]]
double procrustes_inner(
    arma::mat A,
    arma::mat B
) {
    arma::mat M = B * A.t();

    arma::mat u, v;
    arma::vec s;
    arma::svd_econ(u, s, v, M);

    arma::mat R = u * v.t();

    return arma::accu(R % M);
}

//' Bilateral / anisotropic filtering of gradient field
//' 
//' Gradient fields are smoothed using bilateral filtering,
//' in which the smoothed gradient of each edge is computed as
//' the weighted average of the neighboring edges' gradients, considering
//' both distance in space and also similarity in gradients.
//' 
//' @param from_pt,to_pt A pair of `E`-length vectors indicating
//'   the start and end points of each edge (0-indexed).
//' @param field A `2` x `D` x `E` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` latent variables at every edge.
//' @param edges_svd A list with elements `u` (`2` x `2` x `E`), `s` (`2` x `E`), and
//'   `v` (`D` x `2` x `E`) containing the left singular vectors, singular values,
//'   and right singular vectors for each edge.
//' @param coords A `N` x `2` matrix of cell coordinates.
//' @param adj_idx A `N` x `N` sparse adjacency matrix in dgCMatrix format,
//'   containing the mapping from cells to edges. In particular, `adj_idx@x`
//'   should store the 1-indexed edge indices corresponding to each cell-cell pair.
//'   The adjacency matrix is assumed to be symmetric and have zeros on the diagonal.
//'   Importantly, note that the stored edge indices are 1-indexed (*not* 0-indexed here)
//'   in order to avoid problems with R's sparse matrix representation.
//' @param distance Method for computing distance score in weighted average.
//'   See description for details. Defaults to `'euclidean'`.
//' @param similarity Method for computing similarity score in weighted average.
//'   See description for details. Defaults to `'euclidean'`.
//' 
//' @returns A `2` x `D` x `E` array in column-major ordering
//'   containing the smoothed spatial gradient in expression for each of
//'   `D` latent variables at every edge.
//'
// [[Rcpp::export]]
arma::cube smooth_field_edges_cpp(
    arma::uvec & from_pt, // input is 0-indexed
    arma::uvec & to_pt,   // input is 0-indexed
    arma::cube & field,   // arma::cube shape (2, D, E)
    Rcpp::List edges_svd, // u, s, v
    arma::mat & coords,   // arma::mat shape (N, 2)
    const S4& adj_idx,    // dgCMatrix shape (N, N): adj_idx@x should store mapping to 1-indexed edge idx
    unsigned distance,
    unsigned similarity
) {
    // Allocate space for smoothed results
    arma::uword D = field.n_cols;
    arma::uword E = field.n_slices;
    arma::uword N = coords.n_rows;
    arma::cube field_smooth = arma::zeros<arma::cube>(2, D, E);

    // Extract values from adjacency matrix
    if (!Rf_inherits(adj_idx, "dgCMatrix")) throw std::invalid_argument("Adjacency matrix adj_idx must be dgCMatrix");
    arma::uvec adj_i = adj_idx.slot("i");    // row indices (0-based): index of adjacent point
    arma::uvec adj_p = adj_idx.slot("p");    // column pointers (length ncol+1)
    arma::uvec adj_x = adj_idx.slot("x");    // nonzero values: index of corresponding edge
    adj_x -= 1;                              // convert to 0-indexing for edges
    arma::uvec n_neighbors = arma::diff(adj_p);  // number of nonzero values per column

    // Extract singular vectors and values from edge field svd
    arma::cube u = edges_svd["u"];  // (2, 2, E)
    arma::mat s = edges_svd["s"];   // (2, E)
    arma::cube v = edges_svd["v"];  // (D, 2, E)

    for (arma::uword i = 0; i < E; i++) {

        // Edge i: (from -> to)
        arma::uword from = from_pt(i);             // already 0-indexed
        arma::uword to = to_pt(i);                 // already 0-indexed
        arma::rowvec pos_from = coords.row(from);  // length 2 row vector
        arma::rowvec pos_to = coords.row(to);      // length 2 row vector
        arma::mat field_i = field.slice(i);        // 2 x D

        // Allocate weights and values for field weighted average
        arma::vec dvec = arma::zeros<arma::vec>(n_neighbors(from) + n_neighbors(to));           // distance from neighbors
        arma::vec wvec = arma::zeros<arma::vec>(n_neighbors(from) + n_neighbors(to));           // similarity to neighbors
        arma::cube field_mat = arma::zeros<arma::cube>(n_neighbors(from) + n_neighbors(to), 2, D);  // n_neighbors x 2 x D
        arma::rowvec pos_i = arma::zeros<arma::rowvec>(2);

        // Neighboring edge j:
        //   Case 1: (from -> nn_of_from)
        //   Case 2: (nn_of_to -> to)
        arma::uvec nnidx = join_cols(
            adj_i.subvec(adj_p(from), adj_p(from+1)-1),
            adj_i.subvec(adj_p(to), adj_p(to+1)-1)
        );
        arma::uvec neighbor_edge = join_cols(
            adj_x.subvec(adj_p(from), adj_p(from+1)-1),
            adj_x.subvec(adj_p(to), adj_p(to+1)-1)
        );

        for (arma::uword k = 0; k < n_neighbors(from) + n_neighbors(to); k++) {

            // Define unique endpoint of edge i: (from -> to)
            if (k < n_neighbors(from)) {
                pos_i = pos_to;    // Edge j (case 1): (from -> nn_of_from)
            } else {
                pos_i = pos_from;      // Edge j (case 2): (nn_of_to -> to)
            }

            // Extract index and coords of neighboring point of from (case 1) or to (case 2)
            arma::uword nn = nnidx(k);         // neighbor pt idx
            arma::rowvec pos_nn = coords.row(nn);   // length 2 row vector

            // Extract index and field of corresponding edge j
            arma::uword j = neighbor_edge(k);  // edge idx
            arma::mat field_j = field.slice(j);     // 2 x D

            // Compute distance between edges:
            //   Distance is measured between the midpoints of edges
            if (distance == 0) {
                // Euclidean distance:
                //   Absolute distance between the midpoints of edges i and j
                dvec(k) = arma::norm(pos_nn - pos_i) / 2;
            } else if (distance == 1) {
                // Projected distance (anisotropic):
                //   Computes distance between midpoints of i and j projected
                //   along the direction of field_i in coordinate space,
                //   using the left singular vectors of field_i
                dvec(k) = arma::norm(((pos_nn - pos_i) * u.slice(i)).t() % s.col(i)) / 2;
            } else if (distance == 2) {
                dvec(k) = 0;
            } else {
                throw std::invalid_argument("Invalid distance method");
            }
                
            // Compute similarity between fields:
            //   We want to consider only the gradient in the embedding space,
            //   not the directionality in coordinate space, so we solve the
            //   orthogonal Procrustes problem to find the best rotation/reflection
            //   in coordinate space.
            if (similarity == 0) {
                // Euclidean distance:
                //   In the limiting case where field_i and field_j are rank 1,
                //   with embedding gradients z_i and z_j, respectively, then
                //   this is equivalent to computing ||z_i - z_j||. Should equal 0 when i==j.
                arma::mat R = procrustes_mat(field_j, field_i);
                wvec(k) = norm(R * field_j - field_i, "fro");
            } else if (similarity == 1) {
                // Cosine distance (equiv. unit-normalized embedding gradients):
                wvec(k) = std::max(0.0,                      // make non-negative for numerical stability
                    1.0 - procrustes_inner(field_j, field_i)
                        / (arma::norm(s.col(i)) * arma::norm(s.col(j)))
                );
            } else if (similarity == 2) {
                wvec(k) = 0;
            } else {
                throw std::invalid_argument("Invalid similarity method");
            }

            // collect field of neighbor for averaging
            field_mat.row(k) = field_j;
        }
        
        double sig_d = arma::mean(dvec);
        if (sig_d == 0.0) { sig_d = 1; } // avoid division by 0 (if no neighbors or gradient is zero)
        double sig_w = arma::mean(wvec);
        if (sig_w == 0.0) { sig_w = 1; } // avoid division by 0
        arma::vec phi = arma::exp(-(dvec%dvec)/(2*sig_d*sig_d) -(wvec%wvec)/(2*sig_w*sig_w));
        phi = phi / arma::accu(phi);

        for (unsigned d = 0; d < D; d++) {
            field_smooth.subcube(0, d, i, 1, d, i) = field_mat.slice(d).t() * phi;
        }
    }
    
    return field_smooth;
}
