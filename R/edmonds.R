# Vendored from the paper code (BigoneScript.R), lines 463-1165 (marginal cost, Edmonds forest, compare).
# Part of the ZICDT package. Internal engine; user API in zicdt.R.

compute_node_marginal_cost <- function(j, parent_idx, data, penalty_rate, alpha, 
                                       em_max_iter, em_eps) {
  
  # Use C++ function directly - it handles caching internally
  # This ensures consistency between Edmonds and verification
  marginal_score <- compute_node_marginal_cost_cpp(
    j - 1,  # Convert to 0-based indexing
    parent_idx,  # Keep as-is (0 = no parent, 1-based for parent)
    data,
    penalty_rate,
    alpha,
    em_max_iter,
    em_eps
  )
  
  return(marginal_score)
}


# Detect cycle in adjacency matrix and return node indices
# Returns NULL if no cycle, otherwise returns indices of nodes in cycle
detect_cycle_indices <- function(G) {
  n <- nrow(G)
  if (n == 0) return(NULL)
  
  visited <- rep(FALSE, n)
  rec_stack <- rep(FALSE, n)
  
  dfs <- function(v, path = c()) {
    if (rec_stack[v]) {
      # Found cycle - return the cycle node indices
      cycle_start_idx <- which(path == v)
      if (length(cycle_start_idx) > 0) {
        return(path[cycle_start_idx[1]:length(path)])
      }
      return(path)
    }
    if (visited[v]) return(NULL)
    
    visited[v] <<- TRUE
    rec_stack[v] <<- TRUE
    path <- c(path, v)
    
    # Find children of v (vectorized)
    children <- which(G[v, ] == 1)
    for (u in children) {
      cycle <- dfs(u, path)
      if (!is.null(cycle)) return(cycle)
    }
    
    rec_stack[v] <<- FALSE
    return(NULL)
  }
  
  for (v in 1:n) {
    if (!visited[v]) {
      cycle <- dfs(v)
      if (!is.null(cycle)) return(cycle)
    }
  }
  
  return(NULL)
}


# Recursive Edmonds Algorithm Implementation with Proper Cycle Contraction
#
# Correctly handles cycle contraction with proper index tracking and expansion.
# This is the TRUE Edmonds' algorithm for minimum-cost arborescence.
#
# @param cost_matrix Cost matrix (rows = from, cols = to). cost_matrix[i,j] = cost of edge i->j
# @param active_nodes Current active node indices
# @param node_map Mapping from active_nodes to original node indices
# @param original_names Original node names
# @param depth Recursion depth (for logging)
# @return List with G_opt (adjacency) and cost

edmonds_recursive <- function(cost_matrix, active_nodes, node_map, original_names, depth = 0) {
  
  indent <- paste(rep("  ", depth), collapse = "")
  n_active <- length(active_nodes)
  super_root_idx <- active_nodes[1]  # Super-root is first in active_nodes
  real_nodes <- active_nodes[-1]      # All other nodes
  
  cat(sprintf("%sEdmonds iteration (depth %d): %d active nodes\n", 
              indent, depth, n_active - 1))
  
  # BASE CASE: Only super-root remains (all real nodes contracted)
  if (length(real_nodes) == 0) {
    cat(sprintf("%s-> Only super-root remains. Returning empty forest.\n", indent))
    G_final <- matrix(0, nrow(cost_matrix), ncol(cost_matrix))
    return(list(G = G_final, cost = 0))
  }
  
  # ============================================================================
  # PHASE 1: Select minimum incoming edge for each non-root node
  # ============================================================================
  
  min_cost <- rep(Inf, ncol(cost_matrix))
  min_parent <- rep(NA, ncol(cost_matrix))
  
  G_greedy <- matrix(0, nrow(cost_matrix), ncol(cost_matrix))
  
  # Find minimum incoming edge for each real node
  for (j_idx in seq_along(real_nodes)) {
    j <- real_nodes[j_idx]
    
    for (i in active_nodes) {
      if (i != j && is.finite(cost_matrix[i, j]) && cost_matrix[i, j] < min_cost[j]) {
        min_cost[j] <- cost_matrix[i, j]
        min_parent[j] <- i
      }
    }
    
    # Add edge to greedy graph
    if (!is.na(min_parent[j]) && is.finite(min_cost[j])) {
      G_greedy[min_parent[j], j] <- 1
    }
  }
  
  # ============================================================================
  # PHASE 2: Detect cycles (only among real nodes)
  # ============================================================================
  
  G_real <- G_greedy[real_nodes, real_nodes, drop = FALSE]
  cycle_indices <- detect_cycle_indices(G_real)
  
  if (is.null(cycle_indices)) {
    # =========================================================================
    # NO CYCLE: Solution found! Construct final graph
    # =========================================================================
    
    cat(sprintf("%s-> No cycles detected. Solution found!\n", indent))
    
    G_final <- G_greedy
    
    # Compute total cost (sum of all selected node marginal costs)
    total_cost <- 0
    for (j in real_nodes) {
      if (!is.na(min_parent[j]) && is.finite(min_cost[j])) {
        total_cost <- total_cost + min_cost[j]
        if (depth == 0) {
          parent_desc <- if(min_parent[j] == super_root_idx) {
            "SUPER (root)"
          } else {
            paste0("Node ", min_parent[j] - 1)
          }
          cat(sprintf("%s  Node %d: marginal cost = %.6f (parent: %s)\n",
                      indent, j - 1, min_cost[j], parent_desc))
        }
      }
    }
    
    if (depth == 0) {
      cat(sprintf("%s  Edmonds sum of marginal costs: %.8f\n\n", indent, total_cost))
    }
    
    return(list(
      G = G_final,
      cost = total_cost
    ))
  }
  
  # =========================================================================
  # CYCLE FOUND: Apply Edmonds' contraction
  # =========================================================================
  
  # Map cycle indices from G_real (subset) to full active_nodes
  cycle_nodes_active <- real_nodes[cycle_indices]
  
  cycle_names <- original_names[node_map[cycle_nodes_active]]
  cat(sprintf("%s-> Cycle detected: %s\n", indent, paste(cycle_names, collapse = " -> ")))
  
  # =========================================================================
  # PHASE 3: Build contracted graph
  # =========================================================================
  
  # Create mapping for contracted graph
  cycle_representative <- cycle_nodes_active[1]
  nodes_to_remove <- cycle_nodes_active[-1]
  contracted_active <- c(super_root_idx, setdiff(real_nodes, cycle_nodes_active), cycle_representative)
  
  # Map from contracted node index to original node index
  contracted_node_map <- c(
    node_map[super_root_idx],
    node_map[setdiff(real_nodes, cycle_nodes_active)],
    node_map[cycle_representative]
  )
  
  # Build contracted cost matrix
  n_contracted <- length(contracted_active)
  cost_contracted <- matrix(Inf, nrow(cost_matrix), ncol(cost_matrix))
  
  # Copy costs between non-cycle nodes
  for (i in setdiff(contracted_active, cycle_representative)) {
    for (j in setdiff(contracted_active, cycle_representative)) {
      if (i != j) {
        cost_contracted[i, j] <- cost_matrix[i, j]
      }
    }
  }
  
  # Costs FROM external nodes TO cycle (entering the cycle)
  # Formula: cost'(u -> C) = min_v in cycle [ cost(u -> v) - cost(parent(v) -> v) ]
  
  best_entry_node <- rep(NA, length(contracted_active))
  names(best_entry_node) <- as.character(contracted_active)
  
  external_nodes <- setdiff(contracted_active, cycle_representative)
  
  for (u in external_nodes) {
    # Vectorized: find minimum adjusted cost for all cycle nodes
    costs_to_cycle <- cost_matrix[u, cycle_nodes_active]
    parent_costs <- min_cost[cycle_nodes_active]
    adjusted_costs <- costs_to_cycle - parent_costs
    
    valid_mask <- is.finite(adjusted_costs)
    if (any(valid_mask)) {
      min_idx <- which.min(adjusted_costs[valid_mask])
      valid_indices <- which(valid_mask)
      best_v <- cycle_nodes_active[valid_indices[min_idx]]
      min_contracted_cost <- adjusted_costs[valid_mask][min_idx]
      
      cost_contracted[u, cycle_representative] <- min_contracted_cost
      best_entry_node[as.character(u)] <- best_v
    }
  }
  
  # Costs FROM cycle TO external nodes (leaving the cycle)
  # Formula: cost'(C -> w) = min_v in cycle [ cost(v -> w) ]
  
  best_exit_node <- rep(NA, length(contracted_active))
  names(best_exit_node) <- as.character(contracted_active)
  
  for (w in external_nodes) {
    # Vectorized: find minimum cost for all cycle nodes
    costs_from_cycle <- cost_matrix[cycle_nodes_active, w]
    
    valid_mask <- is.finite(costs_from_cycle)
    if (any(valid_mask)) {
      min_idx <- which.min(costs_from_cycle[valid_mask])
      valid_indices <- which(valid_mask)
      best_v <- cycle_nodes_active[valid_indices[min_idx]]
      min_exit_cost <- costs_from_cycle[valid_mask][min_idx]
      
      cost_contracted[cycle_representative, w] <- min_exit_cost
      best_exit_node[as.character(w)] <- best_v
    }
  }
  
  # =========================================================================
  # PHASE 4: Recursively solve on contracted graph
  # =========================================================================
  
  cat(sprintf("%s-> Contracting cycle. New graph size: %d nodes\n", indent, n_contracted - 1))
  
  contracted_names <- original_names[contracted_node_map]
  
  sub_result <- edmonds_recursive(
    cost_contracted,
    contracted_active,
    contracted_node_map,
    contracted_names,
    depth + 1
  )
  
  # =========================================================================
  # PHASE 5: Expand contracted solution back to original graph
  # =========================================================================
  
  cat(sprintf("%s-> Expanding contracted solution\n", indent))
  
  G_final <- matrix(0, nrow(cost_matrix), ncol(cost_matrix))
  
  # First, add all edges from contracted solution that don't involve the cycle
  for (i in setdiff(contracted_active, cycle_representative)) {
    for (j in setdiff(contracted_active, cycle_representative)) {
      if (sub_result$G[i, j] == 1) {
        G_final[i, j] <- 1
      }
    }
  }
  
  # Now handle cycle edges
  # Start with all edges in the original cycle
  for (idx in seq_along(cycle_nodes_active)) {
    v <- cycle_nodes_active[idx]
    parent_v <- min_parent[v]
    if (!is.na(parent_v) && is.finite(min_cost[v])) {
      G_final[parent_v, v] <- 1
    }
  }
  
  # Find incoming edge to cycle in contracted solution
  incoming_to_cycle <- which(sub_result$G[, cycle_representative] == 1)
  
  if (length(incoming_to_cycle) > 0) {
    u_contracted <- incoming_to_cycle[1]
    
    # Find which cycle node this edge targets
    target_v <- best_entry_node[as.character(u_contracted)]
    
    if (!is.na(target_v)) {
      # Remove the original incoming edge to target_v
      parent_target <- min_parent[target_v]
      if (!is.na(parent_target) && is.finite(min_cost[target_v])) {
        G_final[parent_target, target_v] <- 0
      }
      # Add the external edge
      G_final[u_contracted, target_v] <- 1
      cat(sprintf("%s  External edge enters cycle at node %s\n",
                  indent, original_names[node_map[target_v]]))
    }
  } else {
    cat(sprintf("%s  Cycle is independent (all nodes are roots)\n", indent))
  }
  
  # Find outgoing edge from cycle in contracted solution
  outgoing_from_cycle <- which(sub_result$G[cycle_representative, ] == 1)
  
  if (length(outgoing_from_cycle) > 0) {
    for (w_contracted in outgoing_from_cycle) {
      exit_v <- best_exit_node[as.character(w_contracted)]
      if (!is.na(exit_v)) {
        G_final[exit_v, w_contracted] <- 1
        cat(sprintf("%s  Edge leaves cycle from node %s to external node\n",
                    indent, original_names[node_map[exit_v]]))
      }
    }
  }
  
  # =========================================================================
  # PHASE 6: Compute total cost
  # =========================================================================
  
  # Cost accounting:
  # sub_result$cost = cost of contracted solution (with adjusted entry cost)
  # We need to add back:
  #   1. The actual costs of cycle edges
  #   2. The correction for the entry edge (add back what was subtracted)
  
  cat(sprintf("%s-> Cost accounting:\n", indent))
  cat(sprintf("%s  Sub-problem cost: %.6f\n", indent, sub_result$cost))
  
  total_cost <- sub_result$cost
  
  # Add correction for external edge entering the cycle
  incoming_to_cycle <- which(sub_result$G[, cycle_representative] == 1)
  if (length(incoming_to_cycle) > 0) {
    u_contracted <- incoming_to_cycle[1]
    target_v <- best_entry_node[as.character(u_contracted)]
    
    if (!is.na(target_v)) {
      # We used adjusted cost: cost(u->v) - min_cost[v]
      # Now add back min_cost[v] to get the true cost
      correction <- min_cost[target_v]
      cat(sprintf("%s  Entry edge correction: +%.6f (was adjusted by -%.6f)\n",
                  indent, correction, correction))
      total_cost <- total_cost + correction
    }
  }
  
  # Add cost of cycle edges that remain in final solution
  cat(sprintf("%s  Adding cycle edge costs:\n", indent))
  for (v in cycle_nodes_active) {
    parent_v <- min_parent[v]
    if (!is.na(parent_v) && is.finite(min_cost[v]) && G_final[parent_v, v] == 1) {
      cat(sprintf("%s    Cycle edge: node %d -> %d, cost = %.6f\n",
                  indent, parent_v - 1, v - 1, min_cost[v]))
      total_cost <- total_cost + min_cost[v]
    }
  }
  
  cat(sprintf("%s-> Expansion complete. Total cost: %.6f\n", indent, total_cost))
  
  return(list(
    G = G_final,
    cost = total_cost
  ))
}


# Find Optimal Directed Forest Using CORRECT Edmonds' Algorithm
edmonds_forest <- function(data, 
                           penalty_rate = 1.0,
                           alpha = 1.0,
                           em_max_iter = 5000,  # Reduced for faster convergence
                           em_eps = 1e-6,       # Relaxed tolerance
                           precomputed_node_results = NULL) {
  
  n_nodes <- length(data)
  node_names <- names(data)
  if (is.null(node_names)) {
    node_names <- as.character(1:n_nodes)
  }
  
  cat(strrep("=", 80), "\n")
  cat("CORRECT EDMONDS' ALGORITHM FOR DIRECTED FORESTS\n")
  cat("WITH FULL CYCLE CONTRACTION AND PROPER INDEX TRACKING\n")
  cat("GUARANTEED GLOBAL OPTIMUM FOR ARBORESCENCES\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("Number of nodes: %d\n", n_nodes))
  use_precomputed <- !is.null(precomputed_node_results)
  if (use_precomputed) {
    cat("Precomputed node losses provided; skipping EM during marginal cost evaluation.\n")
  } else {
    cat(sprintf("Total marginal cost evaluations: %d\n", 
                n_nodes * (n_nodes + 1)))
  }
  cat(strrep("=", 80), "\n\n")
  if (!use_precomputed) {
    clear_node_cache()  # Clears C++ cache (used for both marginal costs and verification)
  } else {
    if (length(precomputed_node_results) != n_nodes) {
      stop("precomputed_node_results must have one entry per node.")
    }
  }
  
  # Step 1: Compute per-node marginal costs
  cat("Step 1: Computing per-node marginal costs...\n\n")
  
  node_cost_matrix <- matrix(Inf, n_nodes + 1, n_nodes + 1)
  diag(node_cost_matrix) <- Inf
  
  start_time <- proc.time()
  
  for (j in 1:n_nodes) {
    if (use_precomputed) {
      node_entry <- precomputed_node_results[[j]]
      if (is.null(node_entry)) {
        stop(sprintf("Missing precomputed results for node %d", j))
      }
      no_parent_entry <- node_entry[["0"]]
      if (is.null(no_parent_entry)) {
        stop(sprintf("Precomputed results for node %d lack the no-parent case.", j))
      }
      node_cost_matrix[1, j + 1] <- alpha * no_parent_entry$loss
      
      for (i in 1:n_nodes) {
        if (i == j) next
        parent_entry <- node_entry[[as.character(i)]]
        if (!is.null(parent_entry)) {
          node_cost_matrix[i + 1, j + 1] <- alpha * parent_entry$loss + penalty_rate
        }
      }
    } else {
      # Cost of j having no parent
      cost_no_parent <- compute_node_marginal_cost(
        j, 0, data, penalty_rate, alpha, em_max_iter, em_eps
      )
      node_cost_matrix[1, j + 1] <- cost_no_parent
      
      # Cost of j having each possible parent
      for (i in 1:n_nodes) {
        if (i != j) {
          cost_with_parent <- compute_node_marginal_cost(
            j, i, data, penalty_rate, alpha, em_max_iter, em_eps
          )
          node_cost_matrix[i + 1, j + 1] <- cost_with_parent
        }
      }
    }
    cat("\n")
  }
  
  elapsed_time <- (proc.time() - start_time)[3]
  # cat(sprintf("Completed in %.2f seconds\n", elapsed_time))
  # cat(sprintf("C++ cache: %d entries\n\n", get_cache_size()))
  
  # Step 2: Run CORRECT Edmonds' Algorithm
  cat("Step 2: Running Edmonds' Algorithm with full cycle contraction...\n\n")
  
  # Create initial node mapping: each i maps to itself
  initial_node_map <- 1:(n_nodes + 1)
  
  edmonds_result <- edmonds_recursive(
    node_cost_matrix,
    1:(n_nodes + 1),
    initial_node_map,
    c("SUPER", node_names)
  )
  
  G_extended <- edmonds_result$G
  total_cost <- edmonds_result$cost
  
  # Step 3: Extract forest
  cat("\nStep 3: Extracting optimal forest...\n")
  
  G_optimal <- G_extended[2:(n_nodes + 1), 2:(n_nodes + 1)]
  rownames(G_optimal) <- node_names
  colnames(G_optimal) <- node_names
  
  # Find roots
  super_root_children <- which(G_extended[1, 2:(n_nodes + 1)] == 1)
  root_nodes <- super_root_children
  root_names <- node_names[root_nodes]
  
  n_edges <- sum(G_optimal)
  
  cat(sprintf("  Number of edges: %d\n", n_edges))
  cat(sprintf("  Number of roots: %d (%s)\n", 
              length(root_nodes), paste(root_names, collapse = ", ")))
  cat(sprintf("  Structure: %s\n\n", 
              if (length(root_nodes) == 1) "Tree" else "Forest"))
  
  # Step 4: Verify with full DAG evaluation
  cat("Step 4: Verifying total score...\n")
  if (use_precomputed) {
    cat("  Using precomputed node losses (no additional EM evaluations).\n")
    verified_edges <- n_edges
    verified_loss <- 0.0
    for (node in 1:n_nodes) {
      parent_ids <- which(G_optimal[, node] == 1)
      key <- if (length(parent_ids) == 0) "0" else as.character(parent_ids[1])
      entry <- precomputed_node_results[[node]][[key]]
      if (is.null(entry)) {
        stop(sprintf("Missing precomputed results for node %d with parent key %s", node, key))
      }
      verified_loss <- verified_loss + entry$loss
    }
    verified_score <- alpha * verified_loss + penalty_rate * verified_edges
    
    cat(sprintf("  Edmonds total cost: %.8f\n", total_cost))
    cat(sprintf("    = (%.2f x KL loss) + %.2f x %d edges\n", 
                alpha, penalty_rate, n_edges))
    cat(sprintf("  Verified DAG score (precomputed): %.8f\n", verified_score))
    cat(sprintf("    = %.2f x %.8f (KL) + %.2f x %d (edges) = %.8f\n",
                alpha, verified_loss, penalty_rate, verified_edges,
                alpha * verified_loss + penalty_rate * verified_edges))
    
    score_diff <- abs(total_cost - verified_score)
    relative_diff <- score_diff / max(abs(verified_score), 1e-10)
    loss_component_edmonds <- total_cost - penalty_rate * n_edges
    loss_diff <- abs(loss_component_edmonds - verified_loss)
  } else {
    cat(sprintf("Cache before verification: %d entries\n", get_cache_size()))
    
    final_result <- getDAGparametersCpppenalty_cached(
      G_optimal, data, penalty_rate, em_max_iter, em_eps, alpha
    )
    verified_score <- final_result$total_score
    verified_loss <- final_result$total_loss
    verified_edges <- final_result$num_edges
    
    cat(sprintf("Cache hits: %d, misses: %d\n", 
                final_result$cache_hits, final_result$cache_misses))
    
    cat(sprintf("  Edmonds total cost: %.8f\n", total_cost))
    cat(sprintf("    = (%.2f x KL loss) + %.2f x %d edges\n", 
                alpha, penalty_rate, n_edges))
    cat(sprintf("  Verified DAG score: %.8f\n", verified_score))
    cat(sprintf("    = %.2f x %.8f (KL) + %.2f x %d (edges) = %.8f\n",
                alpha, verified_loss, penalty_rate, verified_edges,
                alpha * verified_loss + penalty_rate * verified_edges))
    
    score_diff <- abs(total_cost - verified_score)
    relative_diff <- score_diff / max(abs(verified_score), 1e-10)
    
    # Break down the difference
    loss_component_edmonds <- total_cost - penalty_rate * n_edges
    loss_diff <- abs(loss_component_edmonds - verified_loss)
  }
  
  cat(sprintf("\n  Breakdown:\n"))
  cat(sprintf("    Edmonds KL component: %.8f\n", loss_component_edmonds))
  cat(sprintf("    Verified KL loss:     %.8f\n", alpha * verified_loss))
  cat(sprintf("    KL difference:        %.2e\n", loss_diff))
  
  if (score_diff < 1e-6) {
    cat("\n  [OK] Scores match perfectly!\n")
  } else if (score_diff < 0.01) {
    cat(sprintf("\n  [OK] Scores match within tolerance (diff: %.2e, rel: %.2e)\n", 
                score_diff, relative_diff))
  } else {
    cat(sprintf("\n  ! WARNING: Score mismatch (diff: %.2e, rel: %.2e)\n",
                score_diff, relative_diff))
    if (loss_diff > 1e-6) {
      cat("  -> Issue: KL loss values differ between marginal and full evaluation\n")
      cat("  -> Likely cause: EM convergence or numerical precision\n")
    } else {
      cat("  -> Issue: Penalty calculation differs\n")
      cat("  -> Check edge counting\n")
    }
  }
  
  # Step 5: Sanity checks
  cat("\nStep 5: Sanity checks...\n")
  
  graph_obj <- graph_from_adjacency_matrix(G_optimal, mode = "directed")
  V(graph_obj)$name <- node_names
  
  if (!igraph::is_dag(graph_obj)) {
    stop("ERROR: Result has cycles!")
  } else {
    cat("  [OK] Result is acyclic\n")
  }
  
  in_degrees <- igraph::degree(graph_obj, mode = "in")
  if (any(in_degrees > 1)) {
    stop("ERROR: Some node has >1 parent!")
  } else {
    cat("  [OK] Each node has <=1 parent\n")
  }
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("EDMONDS' ALGORITHM COMPLETED\n")
  cat(strrep("=", 80), "\n")
  
  return(list(
    adjacency_matrix = G_optimal,
    graph = graph_obj,
    total_score = verified_score,
    edmonds_cost = total_cost,
    score_diff = score_diff,
    n_edges = n_edges,
    n_roots = length(root_nodes),
    roots = root_nodes,
    root_names = root_names,
    is_tree = (length(root_nodes) == 1),
    is_forest = TRUE,
    runtime = elapsed_time
  ))
}


# Print cache statistics
print_cache_stats <- function() {
  cat("\nCache Statistics:\n")
  cat(sprintf("  C++ cache: %d entries\n", get_cache_size()))
}


# Visualize Forest Results
visualize_edmonds_forest <- function(result, main = "Edmonds' Optimal Directed Forest") {
  plot(result$graph, 
       main = main,
       edge.arrow.size = 0.5,
       vertex.size = 30,
       vertex.color = "lightblue",
       vertex.label.color = "black",
       edge.color = "darkgray",
       layout = layout_with_sugiyama)
  
  mtext(sprintf("Score: %.2f | Edges: %d | Roots: %s",
                result$total_score,
                result$n_edges,
                paste(result$root_names, collapse = ", ")),
        side = 3, line = 0.5, cex = 0.8)
}


# Compare with ground truth
compare_with_truth <- function(result, true_G) {
  estimated_G <- result$adjacency_matrix
  
  if (nrow(estimated_G) != nrow(true_G) || ncol(estimated_G) != ncol(true_G)) {
    stop("Dimension mismatch between estimated and true graphs")
  }
  
  TP <- sum(estimated_G == 1 & true_G == 1)
  FP <- sum(estimated_G == 1 & true_G == 0)
  FN <- sum(estimated_G == 0 & true_G == 1)
  TN <- sum(estimated_G == 0 & true_G == 0)
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("COMPARISON WITH GROUND TRUTH\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("True Positives:  %d\n", TP))
  cat(sprintf("False Positives: %d\n", FP))
  cat(sprintf("False Negatives: %d\n", FN))
  cat(sprintf("True Negatives:  %d\n", TN))
  cat(sprintf("\nPrecision: %.3f\n", precision))
  cat(sprintf("Recall:    %.3f\n", recall))
  cat(sprintf("F1 Score:  %.3f\n", f1))
  cat(strrep("=", 80), "\n")
  
  return(list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision,
    recall = recall,
    f1 = f1
  ))
}



# Run batch analysis on multiple datasets (PARALLEL VERSION)
# @param data_dir Directory containing dat_*.RData files
# @param dataset_ids Vector of dataset IDs to process (e.g., 1:10)
# @param output_file Path to output text file
# @param penalty_rate Penalty parameter
# @param alpha Alpha parameter
# @param em_max_iter EM max iterations
# @param em_eps EM convergence threshold
# @param create_plots Whether to create and save plots
# @param n_cores Number of cores to use (default: detectCores() - 1)
# @param use_parallel Whether to use parallel processing (default: TRUE)
