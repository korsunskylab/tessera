#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(Rcpp]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Z-score a sparse matrix across each row
//'
//' @param x,p,i Internal data structures from a
//'   sparse matrix in dgCMatrix format.
//' @param ncol,nrow Dimensions of sparse matrix input.
//' @param thresh Z-scores above `thresh` and below
//'   `-thresh` are clipped to `thresh` and `-thresh`,
//'   respectively.
//'
//' @returns A dense matrix in column-major ordering
//'   with dimensions `nrow` x `ncol`.
// [[Rcpp::export]]
arma::mat scaleRows_dgc(const arma::vec& x, const arma::vec& p, const arma::vec& i, 
                        int ncol, int nrow, float thresh) {
    // (0) fill in non-zero elements
    arma::mat res = arma::zeros<arma::mat>(nrow, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res(i[j], c) = x(j);
        }
    }

    // (1) compute means
    arma::vec mean_vec = arma::zeros<arma::vec>(nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            mean_vec(i[j]) += x[j];
        }
    }
    mean_vec /= ncol;
    
    // (2) compute SDs
    arma::vec sd_vec = arma::zeros<arma::vec>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
            nz(i[j])--;
        }
    }
        
    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }
    
    sd_vec = arma::sqrt(sd_vec / (ncol - 1));
    
    // (3) scale values
    res.each_col() -= mean_vec;
    res.each_col() /= sd_vec;
    res.elem(find(res > thresh)).fill(thresh);
    res.elem(find(res < -thresh)).fill(-thresh);
    return res;
}
