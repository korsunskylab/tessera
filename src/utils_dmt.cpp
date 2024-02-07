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

extern unsigned IDX_FROM_PT;
extern unsigned IDX_TO_PT; 
extern unsigned IDX_FROM_TRI; 
extern unsigned IDX_TO_TRI; 
extern unsigned IDX_X0_PT; 
extern unsigned IDX_X1_PT; 
extern unsigned IDX_Y0_PT; 
extern unsigned IDX_Y1_PT; 
extern unsigned IDX_X0_TRI; 
extern unsigned IDX_X1_TRI; 
extern unsigned IDX_Y0_TRI; 
extern unsigned IDX_Y1_TRI; 
extern unsigned IDX_LENGTH_PT; 
extern unsigned IDX_LENGTH_TRI; 
extern unsigned IDX_BOUNDARY; 
extern unsigned IDX_F_PRIM; 
extern unsigned IDX_F_DUAL; 
extern unsigned IDX_AGG_FROM; 
extern unsigned IDX_AGG_TO;


// [[Rcpp::depends(Rcpp]]
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::umat foo_triplets_edges(
    arma::umat & triplets, 
    arma::umat & edges
) {
    unsigned N = edges.n_rows; 
    arma::umat res = arma::zeros<arma::umat>(N, 2); 
    for (unsigned i = 0; i < N; i++) {
        res.row(i) = arma::intersect(triplets.row(edges(i, 0)), triplets.row(edges(i, 1))); 
    }
    return res; 
}
    

// [[Rcpp::export]]
Rcpp::List do_dmt_forest_cpp(
    const arma::vec & f, 
    const arma::uvec & edges_from,
    const arma::uvec & edges_to,
    const arma::vec & edges_f
) { 
    unsigned N = f.n_rows; 
    arma::uvec labels = arma::zeros<arma::uvec>(N); // 0 means unassigned
    std::list<unsigned> critpts; 
    std::list<unsigned> saddles; 
    unsigned forest_size = 0; // number of edges in forest 
    arma::umat forest = arma::zeros<arma::umat>(edges_f.n_rows, 2); // initialize bigger than needed

    // for every point, remember incoming edge 
    // each parent 
    // in a tree, each node has 0 or 1 parent
    arma::uvec parent_edge = arma::zeros<arma::uvec>(N);  
    arma::uvec parent = arma::zeros<arma::uvec>(N);  
    
    unsigned p1, p2, l1, l2, e; 
    arma::uvec edge_order = arma::sort_index(edges_f, "descend"); 
    for (unsigned i = 0; i < edge_order.n_rows; i++) {
        e = edge_order(i); 
        p1 = edges_from(e); 
        p2 = edges_to(e); 
        l1 = labels(p1);
        l2 = labels(p2); 
        
        if (l1 == 0 & l2 == 0) {
            // make new spanning tree
            if (f(p1) > f(p2)) {
                critpts.push_back(p1 + 1); 
                forest(forest_size, 0) = p1; // from 
                forest(forest_size, 1) = p2; // to
                forest_size++; 
                parent_edge(p2) = e; 
                parent(p2) = p1; 
            } else {
                critpts.push_back(p2 + 1); 
                forest(forest_size, 0) = p2; // from 
                forest(forest_size, 1) = p1; // to
                forest_size++; 
                parent_edge(p1) = e; 
                parent(p1) = p2; 
            }
            // labels(p1) = critpts.size(); 
            // labels(p2) = critpts.size(); 
            // label is the index of critical point (1-indexed)
            labels(p1) = critpts.back(); 
            labels(p2) = critpts.back(); 
        } else if (l1 == 0) {
            // continue spanning tree
            labels(p1) = l2; 
            forest(forest_size, 0) = p2; // from 
            forest(forest_size, 1) = p1; // to
            forest_size++; 
            parent_edge(p1) = e; 
            parent(p1) = p2; 
        } else if (l2 == 0) {
            // continue spanning tree
            labels(p2) = l1; 
            forest(forest_size, 0) = p1; // from 
            forest(forest_size, 1) = p2; // to
            forest_size++; 
            parent_edge(p2) = e; 
            parent(p2) = p1; 
        } else if (l1 != l2) {
            // vertices are labeled and belong to different critpts 
            // mark edge as potential saddle
            saddles.push_back(e + 1); 
        } else {
            // do nothing 
        }
    }

    forest = forest.rows(0, forest_size - 1); 
    return Rcpp::List::create(
        Rcpp::Named("edges") = forest + 1, 
        Rcpp::Named("saddles") = saddles, 
        Rcpp::Named("labels") = labels, 
        Rcpp::Named("critpts") = critpts,
        Rcpp::Named("parent") = parent + 1, 
        Rcpp::Named("parent_edge") = parent_edge + 1
    ); 
}


// [[Rcpp::export]]
std::list<unsigned> trace_back_cpp(
    unsigned v0, 
    unsigned vcrit, 
    const arma::uvec & parent_edge, 
    const arma::uvec & parent 
) {
    std::list<unsigned> epath; 
    while (v0 != vcrit) {
        epath.push_back(parent_edge(v0) + 1); 
        v0 = parent(v0);
    }
    return epath; 
}


// [[Rcpp::export]]
std::vector<std::list<unsigned> > trace_epaths_cpp(
    const arma::uvec & saddles, 
    const arma::uvec & vcrits, 
    const arma::uvec & edges_from,  
    const arma::uvec & edges_to, 
    const arma::uvec & parent_edge, 
    const arma::uvec & parent     
) {
    std::vector<std::list<unsigned> > epaths(2 * saddles.n_elem); 
    std::list<unsigned> tmp; 
    for (unsigned i = 0; i < saddles.n_elem; i++) {
        epaths[2 * i] = trace_back_cpp(
            edges_from(saddles(i)), 
            vcrits(edges_from(saddles(i))), 
            parent_edge, 
            parent
        ); 
        epaths[2 * i + 1] = trace_back_cpp(
            edges_to(saddles(i)), 
            vcrits(edges_to(saddles(i))), 
            parent_edge, 
            parent
        ); 
    }
    return epaths; 
}


/* Copy and Pasted from utils_shape.cpp */
arma::uvec arma_setdiff2(const arma::uvec& vec1, const arma::uvec& vec2) {
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

// [[Rcpp::export]]
arma::uvec get_e_sep(
    std::vector<arma::uvec> epaths, // 1-indexed
    arma::uvec saddles, // 1-indexed
    int nedges
) {
    arma::uvec e_sep_bool = arma::zeros<arma::uvec>(nedges); 
    for (int i = 0; i < saddles.n_rows; i++) {
        e_sep_bool(saddles(i) - 1) = 1;         
    }
    for (int i = 0; i < epaths.size(); i++) {
        for (int j = 0; j < epaths[i].n_rows; j++) {
            e_sep_bool(epaths[i](j) - 1) = 1; 
        }
    }
    arma::uvec e_sep = arma::find(e_sep_bool==1); // 0-indexed
    return e_sep; 
}




// [[Rcpp::export]]
arma::uvec prune_e_sep(
    arma::mat & edges, // 1-indexed
    int ntris,
    std::vector<bool>  is_tri_external, // 0-indexed
    arma::uvec e_sep // 0-indexed
) {
    // unsigned IDX_FROM_TRI = 2; 
    // unsigned IDX_TO_TRI = 3; 
    int u, v, e; 
    
    // get all separatrices
    // NOTE: concatenating both vec1 and vec2 into vec right away caused memory issues 
    arma::vec vec1 = edges.submat(e_sep, arma::uvec{IDX_FROM_TRI}) - 1;
    vec1 = arma::unique(vec1); 
    arma::vec vec2 = edges.submat(e_sep, arma::uvec{IDX_TO_TRI}) - 1;
    vec2 = arma::unique(vec2);     
    arma::vec v_sep = arma::zeros<arma::vec>(vec1.n_rows + vec2.n_rows); 
    v_sep.rows(0, vec1.n_rows - 1) = vec1; 
    v_sep.rows(vec1.n_rows, vec1.n_rows + vec2.n_rows - 1) = vec2; 
    v_sep = arma::unique(v_sep); 
    
    // For every tri in v_sep, get its edges on the sep graph 
    std::vector<std::list<int> > tri_to_edge(ntris); 

    int ntot = 0; 
    for (int i = 0; i < e_sep.n_rows; i++) {
        tri_to_edge[edges(e_sep(i), IDX_FROM_TRI) - 1].push_back(e_sep(i)); 
        tri_to_edge[edges(e_sep(i), IDX_TO_TRI) - 1].push_back(e_sep(i)); 
    }

    // stray vertices have degree = 1 and are not on the boundary
    std::stack<int> v_stray; 
    for (int i = 0; i < v_sep.n_rows; i++) {
        v = v_sep[i]; 
        if ((is_tri_external[v]==false) && (tri_to_edge[v].size() == 1)) {
            v_stray.push(v); 
        }
    }

    

    int n_remove = 0; 
    arma::uvec e_sep_remove = arma::zeros<arma::uvec>(e_sep.n_rows); 
    while (!v_stray.empty()) {
       // pop a stray triangle
        v = v_stray.top(); 
        v_stray.pop(); 

        // sometimes, stray nodes are disconnected before we pop them 
        // in this case, just move on to the next one 
        if (tri_to_edge[v].size() == 0) continue;
        
        // v_sep = v_sep(arma::find(v_sep != v)); // not actually needed?
    
        // find its one edge e
        e = tri_to_edge[v].front(); 
        
        // decrease degree of other tri on e 
        if ((edges(e, IDX_FROM_TRI) - 1) == v) {
            u = edges(e, IDX_TO_TRI) - 1; 
        } else {
            u = edges(e, IDX_FROM_TRI) - 1; 
        }
        auto it = std::find(tri_to_edge[u].begin(), tri_to_edge[u].end(), e); 
        tri_to_edge[u].erase(it);
        if (tri_to_edge[u].size() == 1) {
            v_stray.push(u);             
        }
        
        // remove e from separatrix
        e_sep_remove(n_remove++) = e; 
        // e_sep = e_sep(arma::find(e_sep != e)); 
    }
    e_sep = arma_setdiff2(e_sep, e_sep_remove.rows(0, n_remove-1)); 
    return e_sep; 
}

