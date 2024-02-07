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


// arma::vec field_to_f_cpp(
//     arma::cube & field
// ) {
//     unsigned N = field.n_slices; 
//     arma::mat u, v;
//     arma::vec s;
//     arma::vec f = arma::zeros<arma::vec>(N); 
//     for (unsigned i = 0; i < N; i++) {
//         arma::svd(u, s, v, field.slice(i)); 
//         f(i) = arma::accu(s); 
//     }
//     return f; 
// }


// [[Rcpp::export]]
arma::cube smooth_field_cpp(
    arma::uvec & pvec, // diff(adj@p)
    arma::uvec & adj_i, // adj@i
    arma::uvec & adj_p, // adj@p
    arma::cube & field, // arma::cube shape (2, D, N)
    arma::mat & coords, // arma::mat shape (N, 2)
    unsigned distance,
    unsigned similarity
) {
    
    unsigned D = field.n_cols;
    unsigned N = field.n_slices;
    arma::cube field_smooth = arma::zeros<arma::cube>(2, D, N);

    for (unsigned i = 0; i < N; i++) {
        arma::vec pos_a = coords.row(i).t(); // length 2 column vector
        arma::mat field_a = field.slice(i); // 2 x D
        arma::uvec nnidx = adj_i.subvec(adj_p(i), adj_p(i+1)-1); // # neighbor indices (may include self)

        arma::vec dvec = arma::zeros<arma::vec>(pvec(i));    // distance from neighbors
        arma::vec wvec = arma::zeros<arma::vec>(pvec(i));    // similarity to neighbors
        arma::cube field_mat = arma::zeros<arma::cube>(pvec(i), 2, D); // n_neighbors x 2 x D

        for (unsigned j = 0; j < pvec(i); j++) {
            // compute distance between centroids 
            arma::vec pos_b = coords.row(nnidx(j)).t();
            if (distance == 0) {
                // Euclidean distance
                dvec(j) = arma::norm(pos_a - pos_b, 2); 
            } else if (distance == 1) {
                // Projected distance (anisotropic)
                dvec(j) = arma::norm(field_a.t() * (pos_b - pos_a), 2);
            } else if (distance == 2) {
                dvec(j) = 1;
            } else {
                throw std::invalid_argument("Invalid distance method");
            }
                
            // compute similarity between fields
            arma::mat field_b = field.slice(nnidx(j)); // n_coordinates (2) x n_features (D)
            if (similarity == 0) {
                // Euclidean distance
                wvec(j) = arma::norm(field_b - field_a, "fro");
            } else if (similarity == 1) {
                wvec(j) = std::max(0.0, (double) 1 - arma::accu(field_a % field_b) / (arma::norm(field_a, "fro") * arma::norm(field_b, "fro")));
            } else if (similarity == 2) {
                wvec(j) = 1;
            } else {
                throw std::invalid_argument("Invalid similarity method");
            }

            // collect field of neighbor for averaging
            field_mat.row(j) = field_b;
        }
        
        double sig_d = arma::mean(dvec); // problems with division by 0? (if no neighbors)
        if (sig_d == 0.0) { sig_d = 1; } // avoid division by 0
        double sig_w = arma::mean(wvec);
        if (sig_w == 0.0) { sig_w = 1; } // avoid division by 0
        arma::vec phi = arma::exp(-(dvec%dvec)/(2*sig_d*sig_d) -(wvec%wvec)/(2*sig_w*sig_w));
        phi = phi / arma::accu(phi);

        // what if the point has no neighbors?
        for (unsigned k = 0; k < D; k++) {
            field_smooth.subcube(0, k, i, 1, k, i) = field_mat.slice(k).t() * phi;
        }
    }
    
    return field_smooth;
}

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

// [[Rcpp::export]]
arma::cube estimate_field_cpp(
    arma::mat & coords, 
    arma::mat & embeddings, 
    arma::uvec & adj_i, 
    arma::uvec & adj_p 
) {
    unsigned N = coords.n_rows; 
    unsigned D = embeddings.n_cols; 
    arma::cube field = arma::zeros<arma::cube>(2, D, N); 
    // arma::mat field = arma::zeros<arma::mat>(N, 6); 
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


// OLD VERSION 
// arma::mat estimate_field_cpp(
//     arma::mat & coords, 
//     arma::mat & embeddings, 
//     arma::uvec & adj_i, 
//     arma::uvec & adj_p 
    
// ) {
//     unsigned N = coords.n_rows; 
//     arma::mat field = arma::zeros<arma::mat>(N, 6); 
//     arma::uvec nnidx; 
//     arma::mat f, u, v;
//     arma::vec s; 
//     for (unsigned i = 0; i < N; i++) {
//         try {
//             nnidx = adj_i.subvec(adj_p(i), adj_p(i+1)-1); 
//             f = embeddings.rows(nnidx); 
//             arma::svd(u, s, v, arma::cov(coords.rows(nnidx), f)); 
            
//             // since signs are arbitrary, force y>0 for field vectors 
//             if (u(1, 0) < 0) u.col(0) *= -1; 
//             if (u(1, 1) < 0) u.col(1) *= -1; 
            
//             field(i, 0) = s(0); 
//             field(i, 1) = s(1); 
//             field(i, 2) = u(0, 0);
//             field(i, 3) = u(1, 0);
//             field(i, 4) = u(0, 1);
//             field(i, 5) = u(1, 1); 
            
//         } catch (const std::exception& e) {
//             Rcout << "iter " << i << endl;
//             ::Rf_error("error");         
//         }        
        
//     } 
        
//     return field; 
// }



// arma::mat smooth_field_cpp_old(
//     arma::uvec & pvec, // diff(adj@p)
//     arma::uvec & jvec, // adj@i
//     arma::mat & field_coords, 
//     arma::mat & coords,
//     arma::mat & embeddings, 
//     arma::uvec & adj_p, // adj@i
//     bool reorient_grad_pos_y
// ) {
//     arma::mat rotmat = {{0, -1}, {1, 0}}; // check this 
//     unsigned x = 0; 
//     unsigned N = field_coords.n_rows; 
//     arma::mat field_smooth = arma::zeros<arma::mat>(N, 6); 
//     arma::vec pos_a, field_a, pos_b, field_b, dvec, wvec, phi, grad_new, ortho_new; 
//     arma::uvec nnidx; 
//     arma::mat wkmat, coords_nn_scaled; 
//     double sig, len_grad, len_ortho; 
//     for (unsigned i = 0; i < N; i++) {
//         pos_a = coords.row(i).t(); 
//         field_a = field_coords.row(i).t(); 
//         dvec = arma::zeros<arma::vec>(pvec(i)); 
//         wvec = arma::zeros<arma::vec>(pvec(i)); 
//         wkmat = arma::zeros<arma::mat>(pvec(i), 2); 
//         for (unsigned j = 0; j < pvec(i); j++) {
//             // compute distance between centroids 
//             pos_b = coords.row(jvec(x)).t(); 
//             dvec(j) = arma::norm(pos_a - pos_b, 2); 
            
//             // compute dot product between fields
//             field_b = field_coords.row(jvec(x)).t(); 
            
//             wvec(j) = abs(arma::accu(field_a % field_b)); 
//             wkmat.row(j) = (wvec(j) * field_b).t(); // might need to be transpose?
            
//             x++; 
//         }
        
//         sig = arma::mean(dvec); 
//         phi = arma::exp(-(dvec%dvec) / (2 * sig * sig)); 
//         grad_new = (arma::sum(wkmat.each_col() % phi, 0) / arma::accu(wvec % phi)).t(); // dimensions are right here? 
//         grad_new = arma::normalise(grad_new, 2); 
//         ortho_new = rotmat.t() * grad_new; 

//         // keep dy positive for consistency
//         if (reorient_grad_pos_y) {
//             if (grad_new(1) < 0) grad_new *= -1; 
//             if (ortho_new(1) < 0) ortho_new *= -1; 
//         }
    
//         nnidx = jvec.subvec(adj_p(i), adj_p(i+1)-1); // correct indices? 
//         coords_nn_scaled = coords.rows(nnidx); 
//         coords_nn_scaled.each_row() -= arma::mean(coords_nn_scaled, 0); 
//         arma::mat tmp = arma::cov(embeddings.rows(nnidx), coords_nn_scaled * grad_new); 
//         len_grad = sqrt(arma::accu(arma::pow(arma::cov(embeddings.rows(nnidx), coords_nn_scaled * grad_new), 2))); 
//         len_ortho = sqrt(arma::accu(arma::pow(arma::cov(embeddings.rows(nnidx), coords_nn_scaled * ortho_new), 2))); 
    
//         field_smooth(i, 0) = len_grad; 
//         field_smooth(i, 1) = len_ortho; 
//         field_smooth(i, 2) = grad_new(0); 
//         field_smooth(i, 3) = grad_new(1); 
//         field_smooth(i, 4) = ortho_new(0); 
//         field_smooth(i, 5) = ortho_new(1); 
        
//     }
    
//     return field_smooth; 
// }
    

