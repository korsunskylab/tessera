#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <stack>
#include <limits>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(Rcpp]]
// [[Rcpp::depends(RcppArmadillo)]]

unsigned IDX_FROM_PT = 0;
unsigned IDX_TO_PT = 1;
unsigned IDX_FROM_TRI = 2;
unsigned IDX_TO_TRI = 3;
unsigned IDX_X0_PT = 4;
unsigned IDX_X1_PT = 5;
unsigned IDX_Y0_PT = 6;
unsigned IDX_Y1_PT = 7;
unsigned IDX_X0_TRI = 8;
unsigned IDX_X1_TRI = 9;
unsigned IDX_Y0_TRI = 10;
unsigned IDX_Y1_TRI = 11;
unsigned IDX_LENGTH_PT = 12;
unsigned IDX_LENGTH_TRI = 13;
// unsigned IDX_QC = 14;
unsigned IDX_BOUNDARY = 14;
unsigned IDX_F_PRIM = 15;
unsigned IDX_F_DUAL = 16;
unsigned IDX_AGG_FROM = 17;
unsigned IDX_AGG_TO = 18;

arma::uvec arma_setdiff(const arma::uvec& vec1, const arma::uvec& vec2) {
    // Create a temporary vector to store the result
    arma::uvec diff;
    
    // Sort the input vectors to prepare for set difference
    arma::uvec sorted_vec1 = arma::sort(vec1);
    arma::uvec sorted_vec2 = arma::sort(vec2);

    // Indices for iterating over the vectors
    size_t i = 0, j = 0;

    // Iterate through sorted vectors to find the set difference
    while (i < sorted_vec1.n_elem && j < sorted_vec2.n_elem) {
        if (sorted_vec1(i) < sorted_vec2(j)) {
            // Add element to the difference vector if it's in vec1 but not in vec2
            diff.resize(diff.n_elem + 1);
            diff(diff.n_elem - 1) = sorted_vec1(i);
            ++i;
        } else if (sorted_vec1(i) == sorted_vec2(j)) {
            // If the elements are equal, move to the next elements in both vectors
            ++i;
            ++j;
        } else {
            // If the element in vec2 is smaller, move to the next element in vec2
            ++j;
        }
    }

    // Add remaining elements of vec1 (if any) to the difference vector
    while (i < sorted_vec1.n_elem) {
        diff.resize(diff.n_elem + 1);
        diff(diff.n_elem - 1) = sorted_vec1(i);
        ++i;
    }

    return diff;
}

