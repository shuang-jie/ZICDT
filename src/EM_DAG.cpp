// ============================================================================
// EM Algorithm for Compositional Data with DAG Structure Learning
// ============================================================================
// This file contains core functions for:
//   1. EM algorithm for compositional data with parent nodes
//   2. Penalized DAG parameter estimation
//   3. Node-level caching for efficient graph structure search
// ============================================================================

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================================
// EM Algorithm for Compositional Data
// ============================================================================
// Estimates parameters for a compositional target variable given parent nodes
// using the Expectation-Maximization algorithm.
//
// Arguments:
//   x: List of parent node data matrices (each matrix is n x p_j)
//   y: Target node compositional data (n x q matrix)
//   max_iter: Maximum number of EM iterations
//   eps: Convergence tolerance
//
// Returns:
//   List containing:
//     - w: Mixture weights for global + parent components
//     - eta: Global category probabilities
//     - M: List of transition matrices (one per parent)
//     - loglik: Log-likelihood trace
//     - loss: Final KL divergence loss
// ============================================================================
// [[Rcpp::export]]
List EMalgorithm_cpp(const List& x,
                     const arma::mat& y,
                     int max_iter = 5000,  // Reduced from 20000
                     double eps = 1e-6) {  // Relaxed from 1e-8 for faster convergence
    
    // Get dimensions
    const int N = y.n_rows;    // Number of samples
    const int q = y.n_cols;    // Number of categories in target
    const int J = x.length();  // Number of parent nodes
    
    // Store parent matrices and dimensions for faster access
    std::vector<mat> parent_data(J);
    std::vector<int> p(J);  // Dimensions of parent nodes
    for(int j = 0; j < J; j++) {
        parent_data[j] = as<mat>(x[j]);
        p[j] = parent_data[j].n_cols;
    }
    
    // Initialize parameters with data-driven values (better than uniform)
    vec w(J + 1);
    w.fill(1.0/(J + 1));
    
    // Smart initialization: use empirical column means
    vec eta = mean(y, 0).t();  // Much better than uniform 1/q
    
    // Initialize M matrices with smart values (each column = eta for better starting point)
    std::vector<mat> M(J);
    for(int j = 0; j < J; j++) {
        M[j] = mat(q, p[j]);
        for(int c = 0; c < p[j]; c++) {
            M[j].col(c) = eta;  // Initialize each column to empirical marginal
        }
    }
    
    // Pre-allocate memory for E-step (large arrays - use std::vector for performance)
    std::vector<double> gamma0(N * q);      // Global topic responsibility
    std::vector<double> gammaJ(N * q * J);  // Store all parent responsibilities
    std::vector<mat> pi(J);                 // Component probabilities
    
    // Initialize pi matrices
    for(int j = 0; j < J; j++) {
        pi[j] = mat(N * q, p[j]);
    }
    
    // Storage for convergence tracking
    std::vector<double> loglik(max_iter);
    vec cov_contrib(J);  // Small temporary array
    
    // Main EM loop
    double prev_ll = -1e300;  // Initialize to very small number
    double loss = 0.0;
    for(int iter = 0; iter < max_iter; iter++) {
        // Reset matrices - use std::fill for better performance
        std::fill(gamma0.begin(), gamma0.end(), 0.0);
        std::fill(gammaJ.begin(), gammaJ.end(), 0.0);
        
        // Reset log-likelihood for this iteration
        double ll = 0.0;
        loss = 0.0;
        
        // E-Step and log-likelihood calculation
        for(int i = 0; i < N; i++) {
            const int row_offset = i * q;
            for(int r = 0; r < q; r++) {
                // Compute category probability p_ir
                double p_ir = w(0) * eta(r);
                
                // Compute parent contributions
                for(int j = 0; j < J; j++) {
                    cov_contrib(j) = 0.0;
                    const arma::mat& M_j = M[j];
                    const arma::mat& x_j = parent_data[j];
                    
                    for(int c = 0; c < p[j]; c++) {
                        cov_contrib(j) += x_j(i, c) * M_j(r, c);
                    }
                    p_ir += w(j + 1) * cov_contrib(j);
                }
                
                // Add to log-likelihood
                ll += y(i,r) * log(p_ir);
                if(y(i,r) != 0) {
                    loss += y(i,r) * (log(y(i,r)) - log(p_ir));
                }
                
                // Compute responsibilities
                gamma0[row_offset + r] = (w(0) * eta(r)) / p_ir;
                
                // Update parent responsibilities and component probabilities
                for(int j = 0; j < J; j++) {
                    const int j_offset = j * N * q;
                    gammaJ[j_offset + row_offset + r] = (w(j + 1) * cov_contrib(j)) / p_ir;
                    
                    // Update component probabilities
                    mat& pi_j = pi[j];
                    const int base_idx = row_offset + r;
                    
                    if(cov_contrib(j) > eps) {
                        const arma::mat& M_j = M[j];
                        const arma::mat& x_j = parent_data[j];
                        for(int c = 0; c < p[j]; c++) {
                            pi_j(base_idx, c) = (x_j(i, c) * M_j(r, c)) / cov_contrib(j);
                        }
                    } else {
                        const double inv_p = 1.0/p[j];
                        for(int c = 0; c < p[j]; c++) {
                            pi_j(base_idx, c) = inv_p;
                        }
                    }
                }
            }
        }
        
        loglik[iter] = ll;  // Store log-likelihood for this iteration
        
        // Check convergence early (not just after max_iter/2)
        if(iter > 500 && std::abs(ll - prev_ll) < 1e-6) {
            loglik.resize(iter + 1);  // Trim unused space
            break;
        }
        prev_ll = ll;
        
        // M-Step
        // Update w
        w.fill(eps);  // Reset with eps
        
        for(int i = 0; i < N; i++) {
            const int row_offset = i * q;
            for(int r = 0; r < q; r++) {
                w(0) += y(i,r) * gamma0[row_offset + r];
                for(int j = 0; j < J; j++) {
                    const int j_offset = j * N * q;
                    w(j + 1) += y(i,r) * gammaJ[j_offset + row_offset + r];
                }
            }
        }
        
        // Normalize w
        w /= N;
        w /= arma::sum(w);
        
        // Update eta
        eta.fill(eps);  // Reset with eps
        
        for(int r = 0; r < q; r++) {
            for(int i = 0; i < N; i++) {
                eta(r) += y(i,r) * gamma0[i * q + r];
            }
        }
        eta /= arma::sum(eta);
        
        // Update M matrices
        for(int j = 0; j < J; j++) {
            mat& M_j = M[j];
            mat M_new(q, p[j]);
            M_new.fill(eps);
            
            // Compute new values
            for(int r = 0; r < q; r++) {
                for(int c = 0; c < p[j]; c++) {
                    double sum = 0.0;
                    const arma::mat& pi_j = pi[j];
                    const int j_offset = j * N * q;
                    
                    for(int i = 0; i < N; i++) {
                        const int row_offset = i * q;
                        sum += y(i,r) * gammaJ[j_offset + row_offset + r] * pi_j(row_offset + r, c);
                    }
                    M_new(r, c) = sum + eps;
                }
            }
            
            // Column normalize with pre-computed inverse
            for(int c = 0; c < p[j]; c++) {
                double col_sum = 0.0;
                for(int r = 0; r < q; r++) {
                    col_sum += M_new(r, c);
                }
                const double col_sum_inv = 1.0/col_sum;
                for(int r = 0; r < q; r++) {
                    M_j(r, c) = M_new(r, c) * col_sum_inv;
                }
            }
        }
    }
    
    // Convert loglik vector to R object
    NumericVector loglik_r(loglik.begin(), loglik.end());
    
    // Convert M matrices to List for return
    List M_list(J);
    for(int j = 0; j < J; j++) {
        M_list[j] = M[j];
    }
    
    // Return results (w and eta are already arma::vec, can be returned directly)
    return List::create(
        Named("w") = w,
        Named("eta") = eta,
        Named("M") = M_list,
        Named("loglik") = loglik_r,
        Named("loss") = loss
    );
}

// ============================================================================
// DAG Parameter Estimation with Penalty (Non-Cached Version)
// ============================================================================
// Computes parameters and penalized loss for a given DAG structure.
// NOTE: For structure search, use getDAGparametersCpppenalty_cached instead.
//
// Arguments:
//   G: DAG adjacency matrix (n_nodes x n_nodes)
//   y_list: List of compositional data matrices for each node
//   penalty_rate: Penalty coefficient for graph complexity (penalizes edges)
//   em_max_iter: Maximum EM iterations
//   em_eps: EM convergence tolerance
//   alpha: Scaling factor for KL divergence loss
//
// Returns:
//   List containing:
//     - parameters: Node-specific EM results
//     - total_loss: Sum of KL divergence losses
//     - num_edges: Total number of edges in the graph
//     - total_score: alpha * total_loss + penalty_rate * num_edges
// ============================================================================
// [[Rcpp::export]]
List getDAGparametersCpppenalty(const arma::mat& G,
                        const List& y_list,
                        double penalty_rate = 1.0,
                        int em_max_iter = 5000,   
                        double em_eps = 1e-6,      
                        double alpha = 1.0) {      
    
    const int n_nodes = y_list.length();
    List node_parameters(n_nodes);
    double total_loss = 0.0;
    
    // Count edges in graph (sum of all 1's in adjacency matrix)
    int num_edges = accu(G);
    
    // Pre-convert all data matrices once (avoid repeated conversions)
    std::vector<mat> data_matrices(n_nodes);
    for(int i = 0; i < n_nodes; i++) {
        data_matrices[i] = as<mat>(y_list[i]);
    }
    
    // For each node in the DAG
    for(int node = 0; node < n_nodes; node++) {
        const arma::mat& y = data_matrices[node];
        List parent_x;
        std::vector<int> parent_indices;
        
        // Collect parents
        for(int parent = 0; parent < n_nodes; parent++) {
            if(G(parent, node) == 1) {
                parent_x.push_back(data_matrices[parent]);
                parent_indices.push_back(parent);
            }
        }
        
        List node_result;
        if(parent_indices.empty()) {
            // No parents: use colMeans to get eta
            vec eta = mean(y, 0).t();
            
            // Calculate KL divergence loss (vectorized)
            double loss = 0;
            for(uword i = 0; i < y.n_rows; i++) {
                for(uword t = 0; t < y.n_cols; t++) {
                  if(y(i,t) > 0){
                    loss += y(i,t) * (log(y(i,t)) - log(eta(t)));
                  }
                }
            }
            total_loss += loss;

            node_result = List::create(
                Named("type") = "no_parents",
                Named("eta") = eta,
                Named("parents") = parent_indices,
                Named("node_loss") = loss
            );
        } else {
            // Has parents: use EMalgorithm_cpp
            List em_result = EMalgorithm_cpp(parent_x, y, em_max_iter, em_eps);
            double loss = em_result["loss"];
            total_loss += loss;
            
            node_result = List::create(
                Named("type") = "has_parents",
                Named("parents") = parent_indices,
                Named("w") = em_result["w"],
                Named("eta") = em_result["eta"],
                Named("M") = em_result["M"],
                Named("node_loss") = loss
            );
        }
        
        node_parameters[node] = node_result;
    }
    
    // Add node names if they exist in y_list
    if(!Rf_isNull(y_list.names())) {
        node_parameters.names() = y_list.names();
    }
    
    double total_score = alpha * total_loss + penalty_rate * num_edges;

    return List::create(
        Named("G") = G,
        Named("parameters") = node_parameters,
        Named("total_loss") = total_loss,
        Named("num_edges") = num_edges,
        Named("total_score") = total_score
    );
}

// ============================================================================
// Node-Level Caching Infrastructure
// ============================================================================
// Cache stores node_loss for each (node_id, parent_set) pair.
// This dramatically speeds up structure search by avoiding redundant EM fits.
// ============================================================================

// Global cache for node results (persists between function calls)
std::unordered_map<std::string, double> global_node_cache;

// Create unique cache key from node ID and parent set
std::string create_node_cache_key(int node_id, const std::vector<int>& parent_ids) {
    std::string key = "n" + std::to_string(node_id) + "_p";
    std::vector<int> sorted_parents = parent_ids;
    std::sort(sorted_parents.begin(), sorted_parents.end());
    for(int p : sorted_parents) {
        key += "_" + std::to_string(p);
    }
    return key;
}

// ============================================================================
// DAG Parameter Estimation with Caching (Recommended for Structure Search)
// ============================================================================
// Same as getDAGparametersCpppenalty but with node-level caching.
// Stores node_loss for each unique (node, parent_set) configuration.
//
// Usage:
//   - Call clear_node_cache() before starting a new structure search
//   - Use this function during hill-climbing, Edmonds, etc.
//   - Call get_cache_size() to monitor cache growth
//
// Penalty is based on number of edges: penalty_rate * sum(G == 1)
// Cache Key Format: "n{node_id}_p_{parent1}_{parent2}_..." (sorted)
// ============================================================================
// [[Rcpp::export]]
List getDAGparametersCpppenalty_cached(const arma::mat& G,
                                        const List& y_list,
                                        double penalty_rate = 1.0,
                                        int em_max_iter = 5000,
                                        double em_eps = 1e-6,
                                        double alpha = 1.0) {
    
    const int n_nodes = y_list.length();
    double total_loss = 0.0;
    int cache_hits = 0;
    int cache_misses = 0;
    
    // Count edges in graph (sum of all 1's in adjacency matrix)
    int num_edges = accu(G);
    
    // Pre-convert all data matrices once (avoid repeated conversions)
    std::vector<mat> data_matrices(n_nodes);
    for(int i = 0; i < n_nodes; i++) {
        data_matrices[i] = as<mat>(y_list[i]);
    }
    
    // Process each node
    for(int node = 0; node < n_nodes; node++) {
        // Get parents for this node
        std::vector<int> parent_ids;
        for(int parent = 0; parent < n_nodes; parent++) {
            if(G(parent, node) == 1) {
                parent_ids.push_back(parent);
            }
        }
        
        // Create cache key
        std::string cache_key = create_node_cache_key(node, parent_ids);
        
        double node_loss;
        
        // Check cache
        auto cache_it = global_node_cache.find(cache_key);
        if(cache_it != global_node_cache.end()) {
            // Cache hit!
            node_loss = cache_it->second;
            cache_hits++;
        } else {
            // Cache miss - compute
            const arma::mat& y = data_matrices[node];
            
            if(parent_ids.empty()) {
                // No parents: use column means
                vec eta = mean(y, 0).t();
                
                // Calculate KL divergence loss (vectorized)
                node_loss = 0.0;
                for(uword i = 0; i < y.n_rows; i++) {
                    for(uword t = 0; t < y.n_cols; t++) {
                        if(y(i,t) > 0) {
                            node_loss += y(i,t) * (log(y(i,t)) - log(eta(t)));
                        }
                    }
                }
                
            } else {
                // Has parents: use EM
                List parent_x;
                for(int p : parent_ids) {
                    parent_x.push_back(data_matrices[p]);
                }
                
                List em_result = EMalgorithm_cpp(parent_x, y, em_max_iter, em_eps);
                node_loss = as<double>(em_result["loss"]);
            }
            
            // Store in cache
            global_node_cache[cache_key] = node_loss;
            cache_misses++;
        }
        
        total_loss += node_loss;
    }
    
    double total_score = alpha * total_loss + penalty_rate * num_edges;
    
    // Print cache stats every 50 evaluations
    static int eval_count = 0;
    eval_count++;
    if(eval_count % 50 == 0) {
        double hit_rate = cache_hits / (double)(cache_hits + cache_misses) * 100.0;
        Rcout << "Cache: " << cache_hits << " hits, " << cache_misses 
              << " misses (" << hit_rate << "% hit rate, " 
              << global_node_cache.size() << " entries)" << std::endl;
    }
    
    return List::create(
        Named("total_score") = total_score,
        Named("total_loss") = total_loss,
        Named("num_edges") = num_edges,
        Named("cache_hits") = cache_hits,
        Named("cache_misses") = cache_misses
    );
}

// ============================================================================
// COMPUTE KL DIVERGENCE WITH FIXED PARAMETERS (FAST RCPP VERSION)
// ============================================================================

// [[Rcpp::export]]
double compute_kl_with_fixed_params_cpp(const List& test_data_list,
                                       const arma::mat& G_fixed,
                                       const List& learned_params) {
    double total_kl = 0.0;
    
    // Get number of nodes
    const int n_nodes = test_data_list.length();
    
    for (int node = 0; node < n_nodes; node++) {
        mat y = as<mat>(test_data_list[node]);
        const int N = y.n_rows;
        const int q = y.n_cols;
        
        List node_params = learned_params[node];
        std::string type = as<std::string>(node_params["type"]);
        
        if (type == "no_parents") {
            // No parents: use fixed eta
            vec eta = as<vec>(node_params["eta"]);
            
            // Compute KL divergence: sum over samples and categories
            for (int i = 0; i < N; i++) {
                for (int t = 0; t < q; t++) {
                    if (y(i, t) > 0) {
                        // Add small epsilon to eta(t) to avoid log(0) = -Inf
                        // This handles the case where a category never appears in training but appears in test
                        double eta_safe = std::max(eta(t), 1e-10);
                        total_kl += y(i, t) * (log(y(i, t)) - log(eta_safe));
                    }
                }
            }
        } else {
            // Has parents: use fixed EM parameters
            std::vector<int> parent_ids;
            for (int parent = 0; parent < n_nodes; parent++) {
                if (G_fixed(parent, node) == 1) {
                    parent_ids.push_back(parent);
                }
            }
            
            // Get fixed parameters
            vec w = as<vec>(node_params["w"]);
            vec eta = as<vec>(node_params["eta"]);
            List M_list = node_params["M"];
            
            // Convert M matrices to vector of matrices for faster access
            std::vector<mat> M(parent_ids.size());
            for (size_t j = 0; j < parent_ids.size(); j++) {
                M[j] = as<mat>(M_list[j]);
            }
            
            // Compute KL divergence using fixed parameters
            for (int i = 0; i < N; i++) {
                for (int t = 0; t < q; t++) {
                    if (y(i, t) > 0) {
                        // Compute predicted probability using fixed parameters
                        // Follow the same pattern as EMalgorithm_cpp
                        double pred_prob = w(0) * eta(t);
                        
                        // Parent contributions - CORRECT ORDER (sum over c first, then j)
                        for (size_t j = 0; j < parent_ids.size(); j++) {
                            mat parent_y = as<mat>(test_data_list[parent_ids[j]]);
                            const int parent_q = parent_y.n_cols;
                            
                            double cov_contrib = 0.0;
                            for (int c = 0; c < parent_q; c++) {
                                cov_contrib += parent_y(i, c) * M[j](t, c);  // ACCUMULATE first
                            }
                            pred_prob += w(j + 1) * cov_contrib;  // Then use the SUM
                        }
                        
                        // Add KL contribution
                        // Add small epsilon to pred_prob to avoid log(0) = -Inf
                        // This handles the case where predicted probability is zero but test data has positive value
                        double pred_prob_safe = std::max(pred_prob, 1e-10);
                        total_kl += y(i, t) * (log(y(i, t)) - log(pred_prob_safe));
                    }
                }
            }
        }
    }
    
    return total_kl;
}

// ============================================================================
// Marginal Cost Function (For Edmonds Algorithm)
// ============================================================================
// Computes the marginal cost for a single node given a parent (or no parent).
// This ensures Edmonds uses the exact same EM calls as verification.
//
// Returns: alpha * node_loss + penalty_rate * num_edges (0 or 1)
// ============================================================================
// [[Rcpp::export]]
double compute_node_marginal_cost_cpp(int node_idx,
                                      int parent_idx,  // 0 = no parent, >0 = parent index (1-based)
                                      const List& y_list,
                                      double penalty_rate,
                                      double alpha,
                                      int em_max_iter,
                                      double em_eps) {
    
    mat y = as<mat>(y_list[node_idx]);
    
    // Create cache key
    std::vector<int> parent_ids;
    if (parent_idx > 0) {
        parent_ids.push_back(parent_idx - 1);  // Convert to 0-based
    }
    std::string cache_key = create_node_cache_key(node_idx, parent_ids);
    
    double node_loss;
    
    // Check cache
    auto cache_it = global_node_cache.find(cache_key);
    if(cache_it != global_node_cache.end()) {
        node_loss = cache_it->second;
    } else {
        // Cache miss - compute
        if (parent_ids.empty()) {
            // No parent: use column means
            vec eta = mean(y, 0).t();
            
            // Calculate KL divergence loss
            node_loss = 0.0;
            for(uword i = 0; i < y.n_rows; i++) {
                for(uword t = 0; t < y.n_cols; t++) {
                    if(y(i,t) > 0) {
                        node_loss += y(i,t) * (log(y(i,t)) - log(eta(t)));
                    }
                }
            }
        } else {
            // Has parent: use EM
            List parent_x;
            parent_x.push_back(as<mat>(y_list[parent_ids[0]]));
            
            List em_result = EMalgorithm_cpp(parent_x, y, em_max_iter, em_eps);
            node_loss = as<double>(em_result["loss"]);
        }
        
        // Store in cache
        global_node_cache[cache_key] = node_loss;
    }
    
    int num_edges = parent_ids.empty() ? 0 : 1;
    return alpha * node_loss + penalty_rate * num_edges;
}

// ============================================================================
// Cache Management Functions
// ============================================================================

// Clear the global node cache (call before starting a new structure search)
// [[Rcpp::export]]
void clear_node_cache() {
    global_node_cache.clear();
    Rcout << "C++ node cache cleared" << std::endl;
}

// Get the number of cached node configurations
// [[Rcpp::export]]
int get_cache_size() {
    return global_node_cache.size();
}
