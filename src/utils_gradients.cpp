#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(Rcpp]]
// [[Rcpp::depends(RcppArmadillo)]]


// Extract X,Y coordinates from a pts data.frame/data.table (not exported).
static arma::mat get_coords(const List& pts) {
    arma::colvec X = as<arma::colvec>(pts["X"]);
    arma::colvec Y = as<arma::colvec>(pts["Y"]);
    arma::mat coords(X.n_elem, 2);
    coords.col(0) = X;
    coords.col(1) = Y;
    return coords;
}


//' Bilateral / anisotropic filtering of gradient field
//'
//' Gradient fields are smoothed using bilateral filtering,
//' in which the smoothed gradient of each point is computed as
//' the weighted average of the neighbors' gradients, considering
//' both distance in space and also similarity in gradients.
//'
//' @param pts A data.frame or data.table with columns `X` and `Y`
//'   containing the spatial coordinates of each cell.
//' @param pvec,adj_i,adj_p A `N` x `N` sparse adjacency matrix
//'   in dgCMatrix format: `pvec = diff(adj@p)`, `adj_i = adj@i`,
//'   and `adj_p = adj@p`
//' @param field A `2` x `D` x `N` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` latent variables at every point in space.
//' @param distance Method for computing distance score in weighted average.
//'   One of `"euclidean"`, `"projected"`, or `"constant"`.
//' @param similarity Method for computing similarity score in weighted average.
//'   One of `"euclidean"`, `"projected"`, or `"constant"`.
//'
//' @returns A `2` x `D` x `N` array in column-major ordering
//'   containing the smoothed spatial gradient in expression for each of
//'   `D` latent variables at every point in space.
//'
// [[Rcpp::export]]
arma::cube smooth_field_cpp(
    const List pts,         // data.frame/data.table with X,Y columns
    arma::uvec & pvec,      // diff(adj@p)
    arma::uvec & adj_i,     // adj@i
    arma::uvec & adj_p,     // adj@p
    arma::cube & field,     // arma::cube shape (2, D, N)
    std::string distance,
    std::string similarity
) {
    arma::mat coords = get_coords(pts);

    unsigned D = field.n_cols;
    unsigned N = field.n_slices;
    arma::cube field_smooth = arma::zeros<arma::cube>(2, D, N);

    for (unsigned i = 0; i < N; i++) {
        arma::vec pos_a = coords.row(i).t(); // length 2 column vector
        arma::mat field_a = field.slice(i); // 2 x D
        arma::uvec nnidx = adj_i.subvec(adj_p(i), adj_p(i+1)-1); // neighbor indices (may include self)

        arma::vec dvec = arma::zeros<arma::vec>(pvec(i));    // distance from neighbors
        arma::vec wvec = arma::zeros<arma::vec>(pvec(i));    // similarity to neighbors
        arma::cube field_mat = arma::zeros<arma::cube>(pvec(i), 2, D); // n_neighbors x 2 x D

        for (unsigned j = 0; j < pvec(i); j++) {
            // compute distance between centroids
            arma::vec pos_b = coords.row(nnidx(j)).t();
            if (distance == "euclidean") {
                dvec(j) = arma::norm(pos_a - pos_b, 2);
            } else if (distance == "projected") {
                // Projected distance (anisotropic)
                dvec(j) = arma::norm(field_a.t() * (pos_b - pos_a), 2);
            } else if (distance == "constant") {
                dvec(j) = 1;
            } else {
                throw std::invalid_argument("Invalid distance method: " + distance);
            }

            // compute similarity between fields
            arma::mat field_b = field.slice(nnidx(j)); // 2 x D
            if (similarity == "euclidean") {
                wvec(j) = arma::norm(field_b - field_a, "fro");
            } else if (similarity == "projected") {
                wvec(j) = std::max(0.0, (double) 1 - arma::accu(field_a % field_b) / (arma::norm(field_a, "fro") * arma::norm(field_b, "fro")));
            } else if (similarity == "constant") {
                wvec(j) = 1;
            } else {
                throw std::invalid_argument("Invalid similarity method: " + similarity);
            }

            // collect field of neighbor for averaging
            field_mat.row(j) = field_b;
        }

        double sig_d = arma::mean(dvec);
        if (sig_d == 0.0) { sig_d = 1; } // avoid division by 0
        double sig_w = arma::mean(wvec);
        if (sig_w == 0.0) { sig_w = 1; } // avoid division by 0
        arma::vec phi = arma::exp(-(dvec%dvec)/(2*sig_d*sig_d) -(wvec%wvec)/(2*sig_w*sig_w));
        phi = phi / arma::accu(phi);

        for (unsigned k = 0; k < D; k++) {
            field_smooth.subcube(0, k, i, 1, k, i) = field_mat.slice(k).t() * phi;
        }
    }

    return field_smooth;
}

//' Compress a gradient field using SVD
//'
//' Expresses the `2` x `D` total derivative at each location as
//' a pair of `2`-dimensional vectors in the gradient and orthogonal
//' directions.
//'
//' @param field A `2` x `D` x `N` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` latent variables at every point, edge, or triangle.
//'
//' @returns A `N` x `6` matrix with the following attributes for
//'   each location:
//'   \item{dx grad,dy grad}{x,y directions of unit vector in the
//'     direction of greatest change (first singular vector).}
//'   \item{dx ortho,dy ortho}{x,y directions of unit vector orthogonal
//'     to the direction of greatest change (second singular vector).}
//'   \item{|grad|,|ortho|}{Magnitude of directional derivative in the
//'     gradient and orthogonal directions (singular values).}
//'
// [[Rcpp::export]]
arma::mat compress_field_cpp(
    arma::cube & field    
) {
    unsigned N = field.n_slices; 
    arma::mat u, v;
    arma::vec s;
    arma::mat res = arma::zeros<arma::mat>(N, 6); 
    for (unsigned i = 0; i < N; i++) {
        arma::svd(u, s, v, field.slice(i)); 
        res(i, 0) = u(0, 0); // dx grad
        res(i, 1) = u(1, 0); // dy grad
        res(i, 2) = u(0, 1); // dx ortho
        res(i, 3) = u(1, 1); // dy ortho
        res(i, 4) = s(0); // |grad|
        res(i, 5) = s(1); // |ortho|
    }
    return res; 
    
}

//' Compute a spatial gradient field at each point (cell)
//'
//' Distance between neighboring cells is normalized to unit distance
//' so that only the direction from each cell to its neighbors matters.
//' The gradient is then the average gradient in expression of each
//' embedding dimension between the index cell and its neighbors.
//'
//' @param pts A data.frame or data.table with columns `X` and `Y`
//'   containing the spatial coordinates of each cell.
//' @param embeddings A `N` x `D` matrix of cell embeddings.
//' @param adj_i,adj_p A `N` x `N` sparse adjacency matrix
//'   in dgCMatrix format.
//'
//' @returns A `2` x `D` x `N` array in column-major ordering
//'   containing the spatial gradient in expression for each of
//'   `D` embedding dimensions at every point in space.
//'
// [[Rcpp::export]]
arma::cube estimate_field_cpp(
    const List pts,
    arma::mat & embeddings,
    arma::uvec & adj_i,
    arma::uvec & adj_p
) {
    arma::mat coords = get_coords(pts);

    unsigned N = coords.n_rows;
    unsigned D = embeddings.n_cols;
    arma::cube field = arma::zeros<arma::cube>(2, D, N);
    arma::uvec nnidx;
    arma::mat f, xy;
    for (unsigned i = 0; i < N; i++) {
        nnidx = adj_i.subvec(adj_p(i), adj_p(i+1)-1);
        // center the nn PCs
        f = embeddings.rows(nnidx);
        f = f.each_row() - embeddings.row(i);

        // center and scale the XY coordinates
        xy = coords.rows(nnidx);
        xy = xy.each_row() - coords.row(i);
        xy = arma::normalise(xy, 2, 1); // p=2, dim=1 (rows)

        field.slice(i) = (xy.t() * f) / nnidx.n_elem;
    }

    return field;
}

