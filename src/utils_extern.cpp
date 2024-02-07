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

// [[Rcpp::export]]
std::vector<std::list<unsigned> > nn_touches_2_vertices(std::vector<std::list<unsigned> > nn_touches, arma::imat triplets) {
    std::vector<std::list<unsigned> > res(nn_touches.size()); 
    arma::irowvec common; 
    for (unsigned i=0; i < nn_touches.size(); i++) {
        for (const auto& j : nn_touches[i]) {
            common = arma::intersect(triplets.row(i), triplets.row(j-1)); 
            if (common.n_elem == 2) {
                res[i].push_back(j); 
            }            
        }
    }
    return res; 
}

