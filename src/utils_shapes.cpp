#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <stack>
#include <limits>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


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


using namespace Rcpp;
using namespace std;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

// [[Rcpp::depends(Rcpp]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
arma::vec mapToConsecutivePositions(const arma::vec& numbers) {
    // Sort the input numbers
    arma::vec sorted_numbers = arma::sort(numbers);

    // Use an unordered map to count unique numbers and assign consecutive positions
    std::unordered_map<unsigned int, unsigned int> position_map;
    unsigned int count = 1; // Counter for consecutive positions

    // Create a map of numbers to consecutive positions
    for (size_t i = 0; i < sorted_numbers.n_elem; ++i) {
        if (position_map.find(sorted_numbers(i)) == position_map.end()) {
            position_map[sorted_numbers(i)] = count++;
        }
    }

    // Create a vector to store the mapped consecutive positions
    arma::vec mapped_positions(numbers.n_elem);

    // Map the original numbers to their consecutive positions using the map
    for (size_t i = 0; i < numbers.n_elem; ++i) {
        mapped_positions(i) = position_map[numbers(i)];
    }

    return mapped_positions;
}

void reindex_cols(
    arma::mat & edges, 
    int pos1, 
    int pos2 
) {
    arma::vec vals = arma::zeros<arma::vec>(2 * edges.n_rows); 
    vals.rows(0, edges.n_rows-1) = edges.col(pos1); 
    vals.rows(edges.n_rows, 2*edges.n_rows - 1) = edges.col(pos2); 
    vals = mapToConsecutivePositions(vals); 
    edges.col(pos1) = vals.rows(0, edges.n_rows-1); 
    edges.col(pos2) = vals.rows(edges.n_rows, 2*edges.n_rows - 1); 
}        

Graph edges_to_graph(
    arma::mat & edges, 
    int pos1, 
    int pos2
) {
    Graph g;             
    for (int i = 0; i < edges.n_rows; i++) {
        // CAREFUL: edges are 1-indexed
        int source = edges(i, pos1) - 1; 
        int target = edges(i, pos2) - 1; 
        boost::add_edge(source, target, g);
    } 
    return g; 
}


// [[Rcpp::export]]
std::vector<arma::uvec> splitSequence(const arma::vec& y) {
    std::map<double, arma::uvec> groupMap;

    for (unsigned int i = 0; i < y.n_elem; ++i) {
        groupMap[y(i)].insert_rows(groupMap[y(i)].n_rows, arma::uvec{ i });
    }

    std::vector<arma::uvec> groups;
    for (const auto& pair : groupMap) {
        groups.push_back(pair.second);
    }

    return groups;
}

void fleury_algorithm(const Graph& g, arma::uvec& euler_cycle, int start_vertex) {
    Graph g_copy = g;

    euler_cycle = arma::zeros<arma::uvec>(boost::num_vertices(g) * 2); // cycle must have fewer than 2|V| vertices
    int len_cycle = 0; 
    arma::uvec degrees(boost::num_vertices(g_copy), arma::fill::zeros);
    int odd_degree_count = 0;

    for (auto v : boost::make_iterator_range(boost::vertices(g_copy))) {
        degrees[v] = boost::degree(v, g_copy);
        if (degrees[v] % 2 != 0) {
            start_vertex = v;
            odd_degree_count++;
        }
    }

    if (odd_degree_count > 0) {
        Rcpp::stop("No Eulerian cycle exists in the provided graph.");
    }

    std::stack<int> stack;
    int current_vertex = start_vertex;
    stack.push(current_vertex);
    
    while (!stack.empty()) {
        // Rcout << "len_cycle:" << len_cycle << endl;
        // Rcout << "BEGIN" << endl;
        if (degrees[current_vertex] != 0) {
            auto adj_edges = boost::adjacent_vertices(current_vertex, g_copy);
            auto it = adj_edges.first;
            while (it != adj_edges.second) {
                int adjacent_vertex = *it;
                if (degrees[adjacent_vertex] > 0) {
                    break;
                }
                ++it;
            }
            int next_vertex = *it;
            stack.push(next_vertex);
            degrees[current_vertex]--;
            degrees[next_vertex]--;
            boost::remove_edge(current_vertex, next_vertex, g_copy);
            current_vertex = next_vertex;
        } else {
            euler_cycle(len_cycle) = current_vertex; 
            len_cycle++; 
            // euler_cycle = arma::join_cols(euler_cycle, arma::uvec{static_cast<unsigned int>(current_vertex)});
            current_vertex = stack.top();
            stack.pop();
        }
        // euler_cycle.print("cycle:"); 
        // Rcout << "END" << endl;
    }
    // Rcout << "1" << endl;
    euler_cycle = euler_cycle.rows(0, len_cycle-1); 
    // Rcout << "2" << endl;

    // small bug: path has 1st 2 elements repeated
    // put first element at the end 
    euler_cycle = arma::join_cols(euler_cycle.tail(euler_cycle.n_elem - 1), euler_cycle.head(1));     
    // Rcout << "3" << endl;
}


std::vector<int> findLargestComponent(const Graph& g) {
  std::vector<int> component(boost::num_vertices(g));
  int num_components = boost::connected_components(g, &component[0]);

  std::vector<std::vector<int>> components(num_components);
  for (int i = 0; i < component.size(); ++i) {
    components[component[i]].push_back(i);
  }

  auto largest_component = std::max_element(components.begin(), components.end(),
    [](const std::vector<int>& a, const std::vector<int>& b) {
      return a.size() < b.size();
    });

  return *largest_component;
}

// Define a function that initializes a graph from a matrix and finds an Euler cycle
// on the largest component, returning original vertex indices
// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
arma::uvec findLargestComponentEulerCycle(
    arma::mat & edges
) {
  Graph g;

  // Adding edges from the input matrix
  for (int i = 0; i < edges.n_rows; i++) {
    // CAREFUL: edges are 1-indexed
    int source = edges(i, 0) - 1;
    int target = edges(i, 1) - 1; 
    boost::add_edge(source, target, g);
  }

  // Finding the largest component
  std::vector<int> largest_component_indices = findLargestComponent(g);

  arma::uvec euler_cycle;
  // std::vector<int> euler_cycle;
  if (boost::num_vertices(g) == largest_component_indices.size()) {
    // Rcout << "one component: do fleury" << endl;
      // Find Euler cycle in the full graph 
      fleury_algorithm(g, euler_cycle, 0);    
    // Rcout << "one component: do fleury" << endl;
  } else {
    // Rcout << "many components" << endl;
      if (largest_component_indices.empty()) {
        Rcpp::stop("No Eulerian cycle exists in the provided graph.");
      }
    
      // Create a subgraph of the largest component
      Graph largest_graph;
      int v_min = std::numeric_limits<int>::max(); 
      for (auto edge : boost::make_iterator_range(boost::edges(g))) {
        int source = boost::source(edge, g);
        int target = boost::target(edge, g);
        if (
            std::find(largest_component_indices.begin(), largest_component_indices.end(), source) != largest_component_indices.end() &&
            std::find(largest_component_indices.begin(), largest_component_indices.end(), target) != largest_component_indices.end()
        ) {
            if (source < v_min) v_min = source;
            if (target < v_min) v_min = target;
          boost::add_edge(source, target, largest_graph);
        }
      }
    // Rcout << "----1----" << endl;
      // Find Euler cycle in the largest component
      fleury_algorithm(largest_graph, euler_cycle, v_min);    
  }
  // if (euler_cycle.empty()) {
  //   Rcpp::stop("No Eulerian cycle exists in the largest component of the provided graph.");
  // }

    return euler_cycle; 
}


// [[Rcpp::export]]
std::vector<arma::uvec > get_agg_to_boundary_edge(
    arma::mat & E,
    unsigned n_agg
) { 
    
    std::vector<arma::uvec> res(n_agg); 
    for (unsigned i = 0; i < E.n_rows; i++) {    
        if (E(i, IDX_BOUNDARY) == true) { // boundary edge
            res[E(i, IDX_AGG_FROM) - 1].insert_rows(res[E(i, IDX_AGG_FROM) - 1].n_rows, arma::uvec{ i });
            res[E(i, IDX_AGG_TO) - 1].insert_rows(res[E(i, IDX_AGG_TO) - 1].n_rows, arma::uvec{ i });
            // res[E(i, 18) - 1].push_back(i); // agg_from 
            // res[E(i, 19) - 1].push_back(i); // agg_to
        }
    }

    // keep unique elements 
    for (unsigned i = 0; i < n_agg; i++) 
        res[i] = arma::unique(res[i]); 
    
    return res; 
}


// [[Rcpp::export]]
std::vector<arma::uvec> get_agg_to_edge(
    arma::mat & edges, 
    unsigned naggs, 
    bool always_include_boundary
) {
    
    std::vector<arma::uvec> agg_to_edge(naggs); 
    for (unsigned e = 0; e < edges.n_rows; e++) {
        if (edges(e, IDX_AGG_FROM) != edges(e, IDX_AGG_TO)) {
            agg_to_edge[edges(e, IDX_AGG_FROM) - 1].insert_rows(agg_to_edge[edges(e, IDX_AGG_FROM) - 1].n_rows, arma::uvec{ e });
            agg_to_edge[edges(e, IDX_AGG_TO) - 1].insert_rows(agg_to_edge[edges(e, IDX_AGG_TO) - 1].n_rows, arma::uvec{ e });
        } else if (always_include_boundary && edges(e, IDX_BOUNDARY) == true) {
            agg_to_edge[edges(e, IDX_AGG_FROM) - 1].insert_rows(agg_to_edge[edges(e, IDX_AGG_FROM) - 1].n_rows, arma::uvec{ e });        
        }
        
        
        
    }
    return agg_to_edge; 
}


arma::uvec arma_setdiff1(const arma::uvec& vec1, const arma::uvec& vec2) {
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




// CAREFUL: everything is 1-indexed 
// [[Rcpp::export]]
arma::mat get_boundary_graph_cpp(
    arma::uvec e_dual, 
    arma::uvec i_pts, 
    arma::uvec e_prim, 
    arma::mat & E, 
    
    double ntris
) {
    // Rcout << 1 << endl;
    // (0) find edges and points of this aggid 
    arma::uvec e_bridge = arma::intersect(e_dual, e_prim); // 1-indexed
    e_prim = arma_setdiff1(e_prim, e_bridge); 
    unsigned ndual = e_dual.n_elem; 
    unsigned nprim = e_prim.n_elem; 
    unsigned nbridge = e_bridge.n_elem; 
    arma::mat res = arma::zeros<arma::mat>(ndual + nprim + nbridge, 6); 

    // Rcout << 2 << endl;
    // (1) re-index points to avoid collision with triangles
    // dual edges

    
    // arma::uvec i_cols_dual = {2, 3, 8, 10, 9, 11}; // 0-indexed
    arma::uvec i_cols_dual = {IDX_FROM_TRI, IDX_TO_TRI, IDX_X0_TRI, IDX_Y0_TRI, IDX_X1_TRI, IDX_Y1_TRI}; 
    res.rows(0, ndual-1) = E.submat(e_dual, i_cols_dual); 

    // Rcout << 3 << endl;
    // primary edges 
    // e_prim.print("e_prim:");
    // e_bridge.print("e_bridge:");
    // e_dual.print("e_dual:");
    // Rcout << "E dim: " << E.n_rows << "x" << E.n_cols << endl;
    // Rcout << "ntris:" << ntris << endl;
    // arma::uvec i_cols_prim = {0, 1, 4, 6, 5, 7}; // 0-indexed
    arma::uvec i_cols_prim = {IDX_FROM_PT, IDX_TO_PT, IDX_X0_PT, IDX_Y0_PT, IDX_X1_PT, IDX_Y1_PT}; 
    if (nprim > 0) {        
        // some shapes may only have primal edges that are also boundaries
        // i think these singletons 
        res.rows(ndual, ndual+nprim-1) = E.submat(e_prim, i_cols_prim); 
        res.submat(ndual, 0, ndual+nprim-1, 1) += ntris; 
    }

    // Rcout << 4 << endl;
    // (2) create new edges to connect overlapping edges: 
    //     always connects external triangle (to_tri) with one point (to_pt or from_pt)    
    unsigned e, j;
    arma::uvec indices; 
    // Rcout << 5 << endl;
    for (unsigned i = 0; i < nbridge; i++) {
        j = ndual+nprim+i; 
        e = e_bridge(i); 
        res(j, 0) = E(e, IDX_FROM_TRI); 
        res(j, 2) = E(e, IDX_X0_TRI); 
        res(j, 3) = E(e, IDX_Y0_TRI); 

        indices = arma::find(i_pts == (E(e, 0)-1)); // i_pts 0-indexed, E(e, 0) 1-indexed
        if (indices.n_elem > 0) { 
            res(j, 1) = E(e, IDX_FROM_PT) + ntris; 
            res(j, 4) = E(e, IDX_X0_PT); 
            res(j, 5) = E(e, IDX_Y0_PT);    
        } else {
            res(j, 1) = E(e, IDX_TO_PT) + ntris; 
            res(j, 4) = E(e, IDX_X1_PT); 
            res(j, 5) = E(e, IDX_Y1_PT); 
        }
    
        // res(j, 6) = std::sqrt(
        //     (res(j, 4) - res(j, 2)) * (res(j, 4) - res(j, 2)) + 
        //     (res(j, 5) - res(j, 3)) * (res(j, 5) - res(j, 3))
        // ); 
    }

    // Rcout << 6 << endl;
    // map indices to consecutive numbers 
    reindex_cols(res, 0, 1);
    // arma::vec vals = arma::zeros<arma::vec>(2 * res.n_rows); 
    // vals.rows(0, res.n_rows-1) = res.col(0); 
    // vals.rows(res.n_rows, 2*res.n_rows - 1) = res.col(1); 
    // vals = mapToConsecutivePositions(vals); 
    // res.col(0) = vals.rows(0, res.n_rows-1); 
    // res.col(1) = vals.rows(res.n_rows, 2*res.n_rows - 1); 
    
    return res; 

    
}


bool has_euler_cycle(
    Graph & g
) {
    arma::uvec degrees(boost::num_vertices(g), arma::fill::zeros);
    // int odd_degree_count = 0;
    for (auto v : boost::make_iterator_range(boost::vertices(g))) {
        degrees[v] = boost::degree(v, g);
        if (degrees[v] % 2 != 0) {
            return false; 
            
            // odd_degree_count++;
            // if (odd_degree_count > 2) 
                // return false; 
        }
    }
    return true; 
}

// Take a local edge matrix and return xy locations 
arma::mat path_to_xymat(
    arma::mat & edges, // local edges 
    arma::uvec & path, // euler circuit path 
    arma::uvec idx1, // from x0 y0
    arma::uvec idx2 // to x1 y1 
) {
    arma::mat m = arma::zeros<arma::mat>(2 * edges.n_rows, 3); 
    // arma::uvec idx = {0, 4, 6}; // from_pt x0_pt y0_pt
    m.rows(0, edges.n_rows-1) = edges.cols(idx1); 
    // idx = {1, 5, 7}; // to_pt x1_pt y1_pt
    m.rows(edges.n_rows, 2*edges.n_rows-1) = edges.cols(idx2); 
    arma::uvec idx = arma::find_unique(m.col(0)); 
    m = m.rows(idx); 
    idx = sort_index(m.col(0));
    m = m.rows(idx); 
    m = m.rows(path); 
    return m.cols(1, 2); 
}


// [[Rcpp::export]]
Rcpp::List trace_polygons_cpp(
    arma::mat & edges, 
    unsigned naggs, 
    unsigned ntris,
    arma::vec & pts_dmt_component
) {
    // NOTE: edges are 1-indexed
    // Create maps for fast indexing 
    // Rcout << 1 << endl;
    std::vector<arma::uvec> agg_to_pt = splitSequence(pts_dmt_component); // i_pts, 0-indexed
    // Rcout << 2 << endl;
    std::vector<arma::uvec> agg_to_boundary_edge = get_agg_to_boundary_edge(edges, naggs); // e_boundary 0-indexed
    // Rcout << 3 << endl;
    std::vector<arma::uvec> agg_to_edge = get_agg_to_edge(edges, naggs, false); // e_dual 0-indexed

    // Rcout << 4 << endl;
    
    // std::vector<arma::mat> res(naggs); 
    Rcpp::List res(naggs); 
    arma::uvec path; 
    for (int aggid = 0; aggid < naggs; aggid++) {
        // Rcout << aggid << "/" << naggs << endl;

        // if shape has no dual edges
        // then make graph from only primal boundary edges
        if (agg_to_edge[aggid].n_rows == 0) {
            // Case I: Boundary completely in primal graph 
            // Rcout << "I: " << aggid << endl;
            // get boundary edges and graph  
            arma::mat edges_boundary = edges.rows(agg_to_boundary_edge[aggid]);             
            reindex_cols(edges_boundary, IDX_FROM_PT, IDX_TO_PT); // by ref 
            Graph g_boundary = edges_to_graph(edges_boundary, IDX_FROM_PT, IDX_TO_PT);  

            // find Euler circuit 
            fleury_algorithm(g_boundary, path, 0);
            // path = findLargestComponentEulerCycle(edges_boundary); 
            
            res[aggid] = path_to_xymat(
                edges_boundary,
                path, 
                arma::uvec {IDX_FROM_PT, IDX_X0_PT, IDX_Y0_PT}, 
                arma::uvec {IDX_TO_PT, IDX_X1_PT, IDX_Y1_PT}
                
                // arma::uvec {0, 4, 6}, 
                // arma::uvec {1, 5, 7}
            );
            continue; 
        }

        // next, build a graph from dual edges and check if it has an Euler cycle 
        arma::mat edges_dual = edges.rows(agg_to_edge[aggid]); 
        reindex_cols(edges_dual, IDX_FROM_TRI, IDX_TO_TRI); // by ref 
        Graph g_dual = edges_to_graph(edges_dual, IDX_FROM_TRI, IDX_TO_TRI);                 
    
        if (has_euler_cycle(g_dual)) {
            // Case II: Boundary completely in dual graph 
            // Rcout << "II: " << aggid << endl;
            arma::uvec path; 
            fleury_algorithm(g_dual, path, 0);

            res[aggid] = path_to_xymat(
                edges_dual,
                path, 
                arma::uvec {IDX_FROM_TRI, IDX_X0_TRI, IDX_Y0_TRI}, 
                arma::uvec {IDX_TO_TRI, IDX_X1_TRI, IDX_Y1_TRI}
            );   
        } else {
            // Case III: Boundary mix of dual and primal graphs
            // Rcout << "III: " << aggid << endl;
            // agg_to_edge[aggid].print("e_dual:");
            // agg_to_pt[aggid].print("i_pts:");
            // agg_to_boundary_edge[aggid].print("e_prim:");
            arma::mat edges_merged = get_boundary_graph_cpp(
                agg_to_edge[aggid],  // e_dual
                agg_to_pt[aggid], // i_pts
                agg_to_boundary_edge[aggid], // e_prim
                edges, 
                ntris
            ); 
            
            // edges_merged.print("edges_merged:"); 
            
            path = findLargestComponentEulerCycle(edges_merged); 
            // path.print("path:"); 
            res[aggid] = path_to_xymat(
                edges_merged,
                path, 
                arma::uvec {0, 2, 3}, // local indices
                arma::uvec {1, 4, 5} // local indices
            );   
        } 
    // return res; 
        
    }
    return res; 
}


// [[Rcpp::export]]
Rcpp::List trace_polygons_pts_cpp(
    arma::mat & edges, 
    unsigned naggs
) {
    
    /*
        Special case function of the one above. 
        
        Use when dual graph is on points (not triangles) and graph is not pruned at all. 

        Dual graph on unpruned points is by definition closed and never needs to use primary graph. 
    
    */
    
    std::vector<arma::uvec> agg_to_edge = get_agg_to_edge(edges, naggs, true); // e_dual 0-indexed    
    Rcpp::List res(naggs); 
    arma::uvec path; 
    for (int aggid = 0; aggid < naggs; aggid++) {
        // Rcout << aggid << "/" << naggs << endl;
        arma::mat edges_dual = edges.rows(agg_to_edge[aggid]); 
        reindex_cols(edges_dual, IDX_FROM_PT, IDX_TO_PT); // by ref 
        Graph g_dual = edges_to_graph(edges_dual, IDX_FROM_PT, IDX_TO_PT);  

        arma::uvec path; 
        fleury_algorithm(g_dual, path, 0);
        res[aggid] = path_to_xymat(
            edges_dual,
            path, 
            arma::uvec {IDX_FROM_PT, IDX_X0_PT, IDX_Y0_PT}, 
            arma::uvec {IDX_TO_PT, IDX_X1_PT, IDX_Y1_PT}
        );   

    }


    return res; 
}

