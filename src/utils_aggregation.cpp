#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(Rcpp]]
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec st_mats_perimeter(std::vector<arma::mat> & Xvec) {
    arma::vec res = arma::zeros<arma::vec>(Xvec.size());
    arma::mat X; 
    for (int k = 0; k < Xvec.size(); k++) {
        X = Xvec[k]; 
        for (int i = 0; i < X.n_rows-1; i++) {
            res(k) += std::sqrt(
                (X(i, 0) - X(i+1, 0)) * (X(i, 0) - X(i+1, 0)) + 
                (X(i, 1) - X(i+1, 1)) * (X(i, 1) - X(i+1, 1))
            ); 
        }
    }
    return res; 
}


// [[Rcpp::export]]
arma::uvec get_e_minus_epaths_saddles(
    std::vector<arma::uvec> & epaths, // 1-indexed
    arma::uvec saddles, // 1-indexed
    int N
) {
    arma::uvec res = arma::ones<arma::uvec>(N);
    for (int i = 0; i < saddles.n_rows; i++) {
        res(saddles(i) - 1) = 0;         
    }
    for (int i = 0; i < epaths.size(); i++) {
        for (int j = 0; j < epaths[i].n_rows; j++) {
            res(epaths[i](j) - 1) = 0; 
        }
    }
    return arma::find(res==1) + 1; // return 1-indexed
}


// [[Rcpp::export]]
bool is_in_list(
    unsigned target,
    std::list<unsigned> mylist
) {
    return std::find(mylist.begin(), mylist.end(), target) != mylist.end(); 
}


// [[Rcpp::export]]
arma::uvec mergeListsToArmaUVec(const std::list<unsigned>& list1, const std::list<unsigned>& list2) {
    // Calculate the total size needed for the Armadillo uvec
    size_t totalSize = list1.size() + list2.size();

    // Create an Armadillo uvec with the total size
    arma::uvec result(totalSize);

    // Initialize an iterator to copy elements from list1
    auto it = list1.begin();

    // Copy elements from list1 to the Armadillo uvec
    for (size_t i = 0; i < list1.size(); ++i) {
        result(i) = *it;
        ++it;
    }

    // Initialize an iterator to copy elements from list2
    it = list2.begin();

    // Copy elements from list2 to the Armadillo uvec
    for (size_t i = list1.size(); i < totalSize; ++i) {
        result(i) = *it;
        ++it;
    }

    return result;
}

// [[Rcpp::export]]
arma::uvec findDuplicates(const arma::uvec& input) {
    std::unordered_map<unsigned int, int> seen;
    std::vector<unsigned int> duplicates;

    for (const unsigned int& value : input) {
        // Check if the value is already in the map
        if (seen.find(value) != seen.end()) {
            // It's a duplicate, add it to the duplicates vector
            duplicates.push_back(value);
        } else {
            // It's the first occurrence, add it to the map
            seen[value] = 1;
        }
    }

    // Convert the duplicates vector to an Armadillo uvec
    arma::uvec duplicates_uvec = arma::conv_to<arma::uvec>::from(duplicates);

    return duplicates_uvec;
}

// [[Rcpp::export]]
void update_V_cpp(
    arma::mat & V_pcs, 
    arma::vec & V_npts,
    arma::vec & V_perimeter,
    arma::vec & V_area,
    
    
    unsigned e_merge_from, 
    unsigned e_merge_to, 
    double e_merge_edge_length, 
    double e_merge_area, 
    double e_merge_npts, 
    const arma::rowvec & e_merge_pcs, 
    int agg_mode 
) {
    // Update PCs
    double a0, a1, atot;  
    a0 = V_npts(e_merge_from);
    a1 = V_npts(e_merge_to);
    atot = a0 + a1; 
    a0 /= atot; 
    a1 /= atot; 

    unsigned e = e_merge_from; // to save some text 
    // Rcout << "Updating agg=" << e+1 << ", npts from " << V_npts(e) << " to " << e_merge_npts << endl;
    
    // NOTE: update PCA coords but not the LL score because this was already computed for the edge 
    // V_pcs.row(e) = .5 * V_pcs.row(e) + .5 * V_pcs.row(e_merge_to); 
    // V_pcs.row(e) = a0 * V_pcs.row(e) + a1 * V_pcs.row(e_merge_to); 
    V_pcs.row(e) = e_merge_pcs; 
    // V_entropy(e) = e_entropy; 
    
    // V_pcs_norm.row(e_merge_from) = arma::normalise(V_pcs.row(e_merge_from), 2); // p=2 
    // updates by reference 
    V_area(e) = e_merge_area; 
    V_npts(e) = e_merge_npts; 
    
    
    V_perimeter(e) += V_perimeter(e_merge_to) - 2 * e_merge_edge_length; 
    // V_nedges(e) = e_merge_nedges; 
    // V_nedges_internal(e) = e_merge_nedges_internal; 


    switch(agg_mode) {
        case 1: 
            // Static edge scores
            // no need for updates 
            break; 
        case 2: 
            // Simple PCA similarity 
            // do nothing 
            break; 
        case 3: 
            // Merging small aggs with simple PCA similarity 
            // do nothing 
            break; 
        default: 
            Rcpp::stop("Invalid agg_mode"); 
    }

}


// [[Rcpp::export]]
void update_E_cpp(
    arma::mat & V_pcs, 
    arma::vec & V_perimeter,
    arma::vec & V_area,
    arma::vec & V_npts,
    
    
    arma::uvec & E_from,
    arma::uvec & E_to,
    arma::vec & E_npts,
    arma::vec & E_area,
    arma::vec & E_edge_length,
    arma::mat & E_pcs_merge, 
    
    arma::vec & E_w, 
    arma::vec & E_perimeter_merge, 
    arma::vec & E_score_size, 
    
    arma::vec & E_dscore,
    arma::uvec & e_update,

    std::vector<std::list<unsigned> > & V_to_E_from,  
    std::vector<std::list<unsigned> > & V_to_E_to, 
    
    
    double d_mu, 
    double d_sig, 
    int agg_mode,
    double min_npts, 
    double max_npts 
) {

    arma::uvec e_to, e_from;  
    arma::vec d_to, w_to, l_to, d_from, w_from, l_from; 

    unsigned i; 
    arma::vec a0, a1, atot; 
    a0 = V_npts(E_from(e_update)); 
    a1 = V_npts(E_to(e_update)); 
    
    // a0.print("a0:"); 
    // a1.print("a1:"); 
    
    atot = a0 + a1; 
    a0 /= atot; 
    a1 /= atot; 
    

    // update E_pcs
    arma::mat V0_pcs, V1_pcs; 
    V0_pcs = V_pcs.rows(E_from(e_update)); 
    V0_pcs = V0_pcs.each_col() % a0; 
    V1_pcs = V_pcs.rows(E_to(e_update)); 
    V1_pcs = V1_pcs.each_col() % a1; 
    E_pcs_merge.rows(e_update) = V0_pcs + V1_pcs; 
    

    E_perimeter_merge(e_update) = V_perimeter(E_from(e_update)) + V_perimeter(E_to(e_update)) - 2 * E_edge_length(e_update); 

    arma::vec C_from, C_to, C_merge, dC; // handy temp variables for compactness
    switch(agg_mode) {
        case 1: 
            // Static edge weights 
            // do nothing 
            break; 
        case 2: 
            // Simple PCA similarity 
            E_w(e_update) = arma::sqrt(sum(arma::pow(V_pcs.rows(E_from(e_update)) - V_pcs.rows(E_to(e_update)), 2), 1)); 
            E_w(e_update) = .5 -  1 / (1 + arma::exp(-(E_w(e_update) - d_mu) / d_sig)); 
            
            // Penalize scores as they approach size_max 
            // Multiplicative so it doesn't change the sign of E_w
            
            // Rcout << E_score_size(e_update) << endl;
            E_score_size(e_update) = ((-1/max_npts) * V_npts(E_from(e_update)) + 1) % ((-1/max_npts) * V_npts(E_to(e_update)) + 1); 
            
            // Delta compactness penalty
            C_from = a0 % (4 * arma::datum::pi * V_area(E_from(e_update))) / (V_perimeter(E_from(e_update)) % V_perimeter(E_from(e_update))); 
            C_to = a1 % (4 * arma::datum::pi * V_area(E_to(e_update))) / (V_perimeter(E_to(e_update)) % V_perimeter(E_to(e_update))); 
            C_merge = (4 * arma::datum::pi * E_area(e_update)) / (E_perimeter_merge(e_update) % E_perimeter_merge(e_update)); 
            dC = .5 * (C_merge - C_from - C_to + 1); // ranges from 0 to 1

            E_dscore(e_update) = E_w(e_update) % E_score_size(e_update) % dC; 
            
            // set dscore to -Inf if merge agg exceeds size thresholds
            for (i = 0; i < e_update.n_elem; i++) {
                if (E_npts(e_update(i)) >= max_npts) E_dscore(e_update(i)) = -arma::datum::inf; 
            }

            break; 
        case 3:
            /* 
                Merging small outlier aggs
                
                Scoring is the same as case 4
                Instead of enforcing max size, enforce min size 
            */
            
            E_w(e_update) = arma::sqrt(sum(arma::pow(V_pcs.rows(E_from(e_update)) - V_pcs.rows(E_to(e_update)), 2), 1)); 
            // NOTE: only difference from case 4 is "1-" instead of "0.5-"
            //       this makes w always non-negative 
            E_w(e_update) = 1 -  1 / (1 + arma::exp(-(E_w(e_update) - d_mu) / d_sig)); 
            E_score_size(e_update) = ((-1/max_npts) * V_npts(E_from(e_update)) + 1) % ((-1/max_npts) * V_npts(E_to(e_update)) + 1); 
            C_from = a0 % (4 * arma::datum::pi * V_area(E_from(e_update))) / (V_perimeter(E_from(e_update)) % V_perimeter(E_from(e_update))); 
            C_to = a1 % (4 * arma::datum::pi * V_area(E_to(e_update))) / (V_perimeter(E_to(e_update)) % V_perimeter(E_to(e_update))); 
            C_merge = (4 * arma::datum::pi * E_area(e_update)) / (E_perimeter_merge(e_update) % E_perimeter_merge(e_update)); 
            dC = .5 * (C_merge - C_from - C_to + 1); // ranges from 0 to 1

            E_dscore(e_update) = E_w(e_update) % E_score_size(e_update) % dC; 

            // set dscore to -Inf if merge agg exceeds size thresholds
            for (i = 0; i < e_update.n_elem; i++) {
                // Delete this edge if both aggs are already of min size 
                if (V_npts(E_from(e_update(i))) >= min_npts && V_npts(E_to(e_update(i))) >= min_npts) {
                    E_dscore(e_update(i)) = -arma::datum::inf; 
                }
                // Do nothing if at least one agg is too small 
            }
            
            break; 
        default: 
            Rcpp::stop("Invalid agg_mode"); 
    }
    
}




// [[Rcpp::export]]
std::vector<std::list<unsigned> > merge_aggs_cpp(
    arma::mat & V_pcs, 
    arma::vec & V_area,
    arma::vec & V_perimeter,
    arma::vec & V_npts,
    
    
    arma::uvec & E_from, 
    arma::uvec & E_to, 
    arma::vec & E_npts, 
    arma::vec & E_area, 
    arma::vec & E_edge_length, 
    arma::mat & E_pcs_merge, 
    
    arma::vec & E_w, 
    arma::vec & E_perimeter_merge, 
    arma::vec & E_score_size,
    
    arma::vec & E_dscore,
    
    double d_mu, 
    double d_sig, 
    unsigned iter_max,
    int agg_mode,
    double dscore_thresh,
    double min_npts, 
    double max_npts
) {
    // // Load GMM 
    // arma::gmm_diag gmm; 
    // gmm.load(fname_gmm); 
    // int K = gmm.means.n_cols; 
    

    // Instead of merging shapes on the fly, 
    // remember the old IDs of the aggs that were merged into this new agg
    unsigned N = V_npts.n_elem; 
    std::vector<std::list<unsigned> > memory(N); 
    for (unsigned i = 0; i < N; i++) {
        memory[i].push_back(i); 
    }
    
    std::list<unsigned> aggs_remove; 


    // Index: for each vertex, which edges does it belong to
    //        this lets us update edges in place
    //        the fastest optimization step 
    std::vector<std::list<unsigned> > V_to_E_from(N); 
    std::vector<std::list<unsigned> > V_to_E_to(N); 
    for (unsigned e = 0; e < E_from.n_elem; e++) {
        V_to_E_from[E_from(e)].push_back(e); 
        V_to_E_to[E_to(e)].push_back(e); 
    }
    
    unsigned e_max; 
    arma::uvec partners, dup_partners, d_idx, e_update; 
    for (unsigned iter = 0; iter < iter_max; iter++) {        
        // select edge to merge 
        e_max = E_dscore.index_max(); 
        unsigned from = E_from(e_max); 
        unsigned to = E_to(e_max); 
        
        // Rcout << "Removing edge " << e_max + 1 << ": {" << from + 1 << ", " << to + 1 <<  "}" << " (all 1-indexed)" << endl; 
        
        if (E_dscore(e_max) < dscore_thresh) {
            // Rcout << "Breaking on iter " << iter << endl;
            break; 
        }
        E_dscore(e_max) = -arma::datum::inf; 
        aggs_remove.push_back(to); 
    
    
        // Rcout << "start update_V" << endl;
        // (1) Update meta data (V)
        memory[from].merge(memory[to]); 
        update_V_cpp(
            V_pcs, 
            V_npts, 
            V_perimeter, 
            V_area, 
            
            
            from, 
            to, 
            E_edge_length(e_max),
            E_area(e_max),
            E_npts(e_max), 
            E_pcs_merge.row(e_max), 
            agg_mode 
        ); 

        // Rcout << 2 << endl;
    
        // (2) Update edges (E) 
        // reconnect aggs in adjacency graph 
        for (const auto& i : V_to_E_from[to]) 
            E_from(i) = from; 
        V_to_E_from[from].merge(V_to_E_from[to]); 
        V_to_E_from[from].remove(e_max); 
    
        for (const auto& i : V_to_E_to[to]) 
            E_to(i) = from; 
        V_to_E_to[from].merge(V_to_E_to[to]); 
        V_to_E_to[from].remove(e_max); 
    
        // deal with duplicates 
        // partners: vertices connected to .from vertex
        partners = arma::uvec(V_to_E_from[from].size() + V_to_E_to[from].size()); 
        unsigned i = 0; 
        for (const auto& it : V_to_E_from[from]) {
            partners(i) = E_to(it); 
            i++; 
        }
        for (const auto& it : V_to_E_to[from]) {
            partners(i) = E_from(it); 
            i++; 
        }
        dup_partners = findDuplicates(partners); 
        
        // e_update: indices of edges connnected to .from
        //           should be same order as partners 
        // NOTE: V_to_E_ are std::list but I want e_update to be a vector
        e_update = mergeListsToArmaUVec(V_to_E_from[from], V_to_E_to[from]); 
        for (unsigned d = 0; d < dup_partners.n_elem; d++) {
            // CAUTION: i is the local index within partners and dup_partners
            //          not the global index into edges 
            d_idx = arma::find(partners == dup_partners(d)); // always length 2
    
            // Remove duplicate edge i[2] from Index 
            // Remember to need to remove from .from and from partners 
            if (is_in_list(e_update(d_idx(1)), V_to_E_from[from])) {
                V_to_E_from[from].remove(e_update(d_idx(1))); 
                V_to_E_to[dup_partners(d)].remove(e_update(d_idx(1))); 
            } else {
                V_to_E_to[from].remove(e_update(d_idx(1))); 
                V_to_E_from[dup_partners(d)].remove(e_update(d_idx(1))); 
            }
            E_dscore(e_update(d_idx(1))) = -arma::datum::inf; 
            E_edge_length(e_update(d_idx(0))) += E_edge_length(e_update(d_idx(1))); 
            // E_nedges(e_update(d_idx(0))) += E_nedges(e_update(d_idx(1))); 
        }
        // Rcout << 3 << endl;
    
        // same as above, sans duplicates 
        e_update = mergeListsToArmaUVec(V_to_E_from[from], V_to_E_to[from]); 
        
        /* UPDATE EDGE STATS */ 
        for (unsigned e = 0; e < e_update.n_elem; e++) {
            E_npts(e_update(e)) = V_npts(E_from(e_update(e))) + V_npts(E_to(e_update(e)));  
            // E_nedges_internal_merge(e_update(e)) = V_nedges_internal(E_from(e_update(e))) + V_nedges_internal(E_to(e_update(e))) + 1;  
            E_area(e_update(e)) = V_area(E_from(e_update(e))) + V_area(E_to(e_update(e)));  
        }
        // Rcout << "update update_E" << endl;
        /* UPDATE EDGE STATS */ 
        update_E_cpp(
            V_pcs, 
            V_perimeter,
            V_area, 
            V_npts,
            E_from,
            E_to,
            E_npts,
            E_area,
            E_edge_length,
            E_pcs_merge, 
            
            E_w, 
            E_perimeter_merge,             
            E_score_size, 
            
            E_dscore,
            
            
            
            e_update,
            
            V_to_E_from, 
            V_to_E_to, 
            d_mu, 
            d_sig, 
            agg_mode,
            min_npts, 
            max_npts 
        ); 


        
    }
    return memory; 
}







