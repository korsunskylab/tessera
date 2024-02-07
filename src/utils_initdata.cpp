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
arma::uvec assign_unique_rowid_cpp(
    arma::vec X, 
    arma::vec Y
) {
    // first sort by col1 
    arma::uvec o1 = arma::sort_index(X), o2;
    int val = X(o1(0)); 
    unsigned loc1; 
    loc1 = 0; 
    for (int i = 0; i < o1.n_rows; i++) {
        if (val != X(o1(i))) {
            // we hit a new value of X
            // go back and get order of these values by Y
            o2 = o1.rows(loc1, i-1); 
            o2 = o2.rows(arma::sort_index(Y(o1.rows(loc1, i-1)))); 
            o1.rows(loc1, i-1) = o2; 
            loc1 = i; 
            val = X(o1(i)); 
        }
    }
    
    // after you leave the loop, remember to handle the last range 
    int i = o1.n_rows; 
    o2 = o1.rows(loc1, i-1); 
    o2 = o2.rows(arma::sort_index(Y(o1.rows(loc1, i-1)))); 
    o1.rows(loc1, i-1) = o2; 
    
    
    int count = 0; 
    // then iterate until you hit a new value
    int val1 = X(o1(0)), val2 = Y(o1(0)); 
    arma::uvec res = arma::zeros<arma::uvec>(o1.n_rows); 
    for (int i = 0; i < o1.n_rows; i++) {
        if (X(o1(i)) == val1) {
            if (Y(o1(i)) == val2) {
                res(o1(i)) = count; 
            } else {
                res(o1(i)) = ++count; 
                val2 = Y(o1(i)); 
            }
        } else {
            res(o1(i)) = ++count; 
            val1 = X(o1(i)); 
            val2 = Y(o1(i)); 
        }
    }
    return res; 
    
}

// [[Rcpp::export]]
arma::mat init_edges_cpp(
    arma::mat & triplets, 
    arma::mat & pts, 
    arma::mat & tris
) {

    // First, take all triplets and make edges from them 
    arma::mat edges = arma::zeros<arma::mat>(triplets.n_rows * 3, 14); 
    arma::rowvec x; 
    for (int i = 0; i < triplets.n_rows; i++) {
        x = arma::sort(triplets.row(i), "ascend"); // preserve from < to 
        edges((i*3), 0) = x(0); 
        edges((i*3), 1) = x(1); 
        edges((i*3) + 1, 0) = x(0); 
        edges((i*3) + 1, 1) = x(2); 
        edges((i*3) + 2, 0) = x(1); 
        edges((i*3) + 2, 1) = x(2); 
    }

    // These triplets contain make non-unique combinations 
    // Use hash function to find duplicates 
    // Also make hash sequential (0-indexed)
    arma::uvec hash = assign_unique_rowid_cpp(edges.col(0), edges.col(1));  // already 0-indexed

    // only keep unique hashes 
    // order edges by their hash values
    arma::uvec i_unique = arma::find_unique(hash); 
    arma::uvec o = arma::sort_index(hash(i_unique)); 
    edges = edges.rows(i_unique); 
    // make edges ordered by hash value
    edges = edges.rows(o); 

    
    // triplets_e: triplet edges 
    // each triplet should be the edges that belong to one triangle
    // does not need to be in same order as triplets input
    arma::mat triplets_e(hash.n_rows / 3, 3); 
    for (int i = 0; i < hash.n_rows; i++) {
        triplets_e(std::floor(i / 3), i % 3) = hash(i); 
    }
    
    int u; 
    for (int e = 0; e < triplets_e.n_rows; e++) {
        for (int i = 0; i < 3; i++) {
            u = triplets_e(e, i); 
            if (edges(u, 3) == 0) {
                // to is empty, put it here 
                edges(u, 3) = e + 1; 
            } else if (edges(u, 2) == 0) {
                // to is full, from is empty 
                if (e + 1 < edges(u, 3)) {
                     // to is correct, fill from
                    edges(u, 2) = e + 1; 
                } else {
                     // swap from and to 
                     // fill from 
                    edges(u, 2) = edges(u, 3);
                    edges(u, 3) = e + 1; 
                }
            } 
        }
    }
    
    for (int e = 0; e < edges.n_rows; e++) {
        if (edges(e, 2) == 0) // from_tri
            edges(e, 2) = arma::datum::nan; 
    }


    edges.col(4) = pts.submat(arma::conv_to<arma::uvec>::from(edges.col(0)-1), arma::uvec {0}); // x0_pt
    edges.col(5) = pts.submat(arma::conv_to<arma::uvec>::from(edges.col(1)-1), arma::uvec {0}); // x1_pt
    edges.col(6) = pts.submat(arma::conv_to<arma::uvec>::from(edges.col(0)-1), arma::uvec {1}); // y0_pt
    edges.col(7) = pts.submat(arma::conv_to<arma::uvec>::from(edges.col(1)-1), arma::uvec {1}); // y1_pt
    
    edges.col(8) = tris.submat(arma::conv_to<arma::uvec>::from(edges.col(2)-1), arma::uvec {0}); // x0_tri
    edges.col(9) = tris.submat(arma::conv_to<arma::uvec>::from(edges.col(3)-1), arma::uvec {0}); // x1_tri
    edges.col(10) = tris.submat(arma::conv_to<arma::uvec>::from(edges.col(2)-1), arma::uvec {1}); // y0_tri
    edges.col(11) = tris.submat(arma::conv_to<arma::uvec>::from(edges.col(3)-1), arma::uvec {1}); // y1_tri

    arma::uvec i_na = arma::find_nan(edges.col(2));
    for (int i = 0; i < i_na.n_rows; i++) {
        edges(i_na(i), 8) = arma::datum::nan; 
        edges(i_na(i), 10) = arma::datum::nan; 
    }
    
    // length_pt
    edges.col(12) = arma::sqrt(
        (edges.col(5) - edges.col(4)) % (edges.col(5) - edges.col(4)) + 
        (edges.col(7) - edges.col(6)) % (edges.col(7) - edges.col(6))
    );
    
    // length_tri
    edges.col(13) = arma::sqrt(
        (edges.col(9) - edges.col(8)) % (edges.col(9) - edges.col(8)) + 
        (edges.col(11) - edges.col(10)) % (edges.col(11) - edges.col(10))
    );
    
    return edges; 
}
    


// [[Rcpp::export]]
arma::mat init_tris_cpp(
    arma::umat & triplets, 
    arma::mat & pts
) {
    arma::mat tris = arma::zeros<arma::mat>(triplets.n_rows, 3); 
    // mean X 
    tris.col(0) = (
        pts.submat(triplets.col(0)-1, arma::uvec {0}) + 
        pts.submat(triplets.col(1)-1, arma::uvec {0}) + 
        pts.submat(triplets.col(2)-1, arma::uvec {0})         
    );
    // mean Y
    tris.col(1) = (
        pts.submat(triplets.col(0)-1, arma::uvec {1}) + 
        pts.submat(triplets.col(1)-1, arma::uvec {1}) + 
        pts.submat(triplets.col(2)-1, arma::uvec {1})         
    );
    tris /= 3; 

    // Area using Herons formula 
    double a, b, c, s; 
    for (int i = 0; i < triplets.n_rows; i++) {
        a = std::sqrt(
            std::pow(pts(triplets(i, 0)-1, 0) - pts(triplets(i, 1)-1, 0), 2) + 
            std::pow(pts(triplets(i, 0)-1, 1) - pts(triplets(i, 1)-1, 1), 2)  
        );
        b = std::sqrt(
            std::pow(pts(triplets(i, 0)-1, 0) - pts(triplets(i, 2)-1, 0), 2) + 
            std::pow(pts(triplets(i, 0)-1, 1) - pts(triplets(i, 2)-1, 1), 2)    
        );
        c = std::sqrt(
            std::pow(pts(triplets(i, 1)-1, 0) - pts(triplets(i, 2)-1, 0), 2) + 
            std::pow(pts(triplets(i, 1)-1, 1) - pts(triplets(i, 2)-1, 1), 2)    
        );
        
        s = .5 * (a + b + c);
        tris(i, 2) = std::sqrt(s * (s-a) * (s-b) * (s-c));
    }
    return tris; 
}


    
