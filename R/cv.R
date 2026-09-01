# Vendored from the paper code (BigoneScript.R), lines 1491-2028 (cross-validation + alpha selection).
# Part of the ZICDT package. Internal engine; user API in zicdt.R.

create_cv_folds <- function(data_list, K = 10) {
  N <- nrow(data_list[[1]])  # Number of samples (same for all nodes)
  fold_ids <- sample(rep(1:K, length.out = N))
  return(fold_ids)
}

# Split list-based data for cross-validation
# @param data_list List of p matrices, each N x q_j
# @param fold_ids Vector of fold assignments
# @param test_fold Which fold to use as test set
# @return List with train and test data
split_data_cv <- function(data_list, fold_ids, test_fold) {
  train_mask <- fold_ids != test_fold
  test_mask <- fold_ids == test_fold
  
  # Split each node's data matrix
  train_data <- lapply(data_list, function(node_data) {
    node_data[train_mask, , drop = FALSE]
  })
  
  test_data <- lapply(data_list, function(node_data) {
    node_data[test_mask, , drop = FALSE]
  })
  
  return(list(train = train_data, test = test_data))
}

# Precompute per-node losses and EM parameters for all possible parents
# @param data_list List of node data matrices (training split)
# @param em_max_iter Maximum EM iterations
# @param em_eps Convergence tolerance for EM
# @param allow_parallel Whether to parallelize over nodes (Unix-like systems only)
# @param n_cores Number of cores for node-level parallelism (NULL to auto-select)
# @return List indexed by node where each entry contains losses/parameters for parent choices
precompute_node_results <- function(data_list,
                                    em_max_iter = 20000,
                                    em_eps = 1e-8,
                                    allow_parallel = FALSE,
                                    n_cores = NULL) {
  n_nodes <- length(data_list)
  node_indices <- seq_len(n_nodes)
  
  compute_for_node <- function(node_idx) {
    node_matrix <- as.matrix(data_list[[node_idx]])
    node_results <- list()
    
    # No parent case (super-root)
    eta <- colMeans(node_matrix)
    loss_no_parent <- 0.0
    for (i in seq_len(nrow(node_matrix))) {
      row_vals <- node_matrix[i, ]
      for (t in seq_len(ncol(node_matrix))) {
        val <- row_vals[t]
        if (val > 0) {
          loss_no_parent <- loss_no_parent + val * (log(val) - log(eta[t]))
        }
      }
    }
    node_results[["0"]] <- list(
      loss = loss_no_parent,
      params = list(
        type = "no_parents",
        parent = 0,
        eta = eta
      )
    )
    
    # Single-parent cases
    for (parent_idx in node_indices) {
      if (parent_idx == node_idx) next
      parent_matrix <- as.matrix(data_list[[parent_idx]])
      em_result <- EMalgorithm_cpp(list(parent_matrix), node_matrix, em_max_iter, em_eps)
      loss_parent <- as.numeric(em_result$loss)
      node_results[[as.character(parent_idx)]] <- list(
        loss = loss_parent,
        em_iters = length(em_result$loglik),   # actual EM iterations run (== em_max_iter means hit cap)
        params = list(
          type = "with_parents",
          parent = parent_idx,
          w = em_result$w,
          eta = em_result$eta,
          M = em_result$M
        )
      )
    }
    
    node_results
  }
  
  use_parallel <- allow_parallel && n_nodes > 1 && .Platform$OS.type != "windows"
  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- min(parallel::detectCores(), n_nodes)
    }
    precomputed <- parallel::mclapply(node_indices, compute_for_node, mc.cores = n_cores)
  } else {
    precomputed <- lapply(node_indices, compute_for_node)
  }
  
  names(precomputed) <- if (!is.null(names(data_list))) names(data_list) else as.character(node_indices)
  precomputed
}

# Train model with alpha penalty and save learned parameters
# @param train_data_list List of training data matrices
# @param alpha Alpha hyperparameter value for training
# @param precomputed_results Optional precomputed node results to reuse EM outputs
# @param penalty_rate Penalty parameter (default 1.0)
# @param em_max_iter EM maximum iterations (default 20000)
# @param em_eps EM convergence tolerance (default 1e-8)
# @return List with model result, learned parameters, and (optionally) reused precomputed data
train_model_cv <- function(train_data_list,
                           alpha,
                           precomputed_results = NULL,
                           penalty_rate = 1.0,
                           em_max_iter = 20000,
                           em_eps = 1e-8) {
  if (is.null(precomputed_results)) {
    precomputed_results <- precompute_node_results(
      train_data_list,
      em_max_iter = em_max_iter,
      em_eps = em_eps
    )
  }
  
  # Train WITH alpha penalty (consistent with evaluation)
  result <- edmonds_forest(train_data_list, 
                           penalty_rate = penalty_rate,
                           alpha = alpha,
                           em_max_iter = em_max_iter,
                           em_eps = em_eps,
                           precomputed_node_results = precomputed_results)
  
  # Extract learned parameters for each node directly from precomputed results
  G_learned <- result$adjacency_matrix
  learned_params <- vector("list", length = nrow(G_learned))
  
  for (node in seq_len(nrow(G_learned))) {
    parent_ids <- which(G_learned[, node] == 1)
    if (length(parent_ids) == 0) {
      entry <- precomputed_results[[node]][["0"]]
    } else {
      parent_id <- parent_ids[1]
      entry <- precomputed_results[[node]][[as.character(parent_id)]]
      if (is.null(entry)) {
        stop(sprintf("Precomputed parameters missing for node %d with parent %d.", node, parent_id))
      }
    }
    learned_params[[node]] <- entry$params
  }
  
  if (!is.null(names(train_data_list))) {
    names(learned_params) <- names(train_data_list)
  }
  
  return(list(
    model = result,
    learned_params = learned_params,
    precomputed = precomputed_results
  ))
}

# Evaluate model on test data using fixed training parameters (true CV)
# @param training_result Result from train_model_cv (contains model + learned_params)
# @param test_data_list List of test data matrices
# @param alpha Alpha hyperparameter (not used in evaluation, just for compatibility)
# @return KL loss using fixed training parameters (true generalization measure)
evaluate_model_cv <- function(training_result, test_data_list, alpha) {
  # Get learned structure and parameters from training
  G_learned <- training_result$model$adjacency_matrix
  learned_params <- training_result$learned_params
  
  # Compute KL loss using FIXED training parameters (no re-fitting on test data)
  # Use fast C++ implementation
  test_kl_loss <- compute_kl_with_fixed_params_cpp(test_data_list, G_learned, learned_params)
  
  return(test_kl_loss)
}

# Cross-validate alpha selection
# @param data_list List of p matrices, each N x q_j
# @param alpha_grid Vector of alpha values to test
# @param K Number of CV folds (default: 5)
# @param use_parallel Whether to parallelize across folds (default: TRUE)
# @param n_cores Number of cores for fold-level parallelism (default: detectCores() - 1)
# @param precompute_parallel Whether to parallelize node precomputation within each fold (Unix-like only)
# @param precompute_n_cores Number of cores for node-level precomputation (NULL to auto-select)
# @param em_max_iter EM maximum iterations for training/precomputation
# @param em_eps EM convergence tolerance for training/precomputation
# @return Data frame with CV results
cross_validate_alpha <- function(data_list,
                                 alpha_grid,
                                 K = 10,
                                 use_parallel = TRUE,
                                 n_cores = NULL,
                                 precompute_parallel = FALSE,
                                 precompute_n_cores = NULL,
                                 em_max_iter = 20000,
                                 em_eps = 1e-8) {
  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    n_cores <- max(1, min(n_cores, K))
  } else {
    n_cores <- 1
  }
  
  if (precompute_parallel && .Platform$OS.type == "windows") {
    warning("Node-level precomputation parallelism is not supported on Windows; falling back to sequential precomputation.")
    precompute_parallel <- FALSE
  }
  
  cat("================================================================================\n")
  cat("CROSS-VALIDATION FOR ALPHA SELECTION\n")
  cat("================================================================================\n")
  cat(sprintf("Data: %d samples, %d nodes\n", nrow(data_list[[1]]), length(data_list)))
  cat(sprintf("Alpha grid: [%s]\n", paste(alpha_grid, collapse = ", ")))
  cat(sprintf("CV folds: %d\n", K))
  cat(sprintf("Parallel across folds: %s\n", if (use_parallel && n_cores > 1) sprintf("YES (%d cores)", n_cores) else "NO"))
  cat(sprintf("Node precomputation parallelism: %s\n", if (precompute_parallel) "YES" else "NO"))
  cat("Precomputing node losses once per fold and reusing them across all alpha values.\n")
  cat("================================================================================\n")
  
  fold_ids <- create_cv_folds(data_list, K)
  folds <- seq_len(K)
  
  process_single_fold <- function(fold) {
    cat(sprintf("-> Starting fold %d of %d\n", fold, K), file = stderr())
    data_split <- split_data_cv(data_list, fold_ids, fold)
    
    train_precomputed <- precompute_node_results(
      data_split$train,
      em_max_iter = em_max_iter,
      em_eps = em_eps,
      allow_parallel = precompute_parallel,
      n_cores = precompute_n_cores
    )
    
    fold_rows <- vector("list", length(alpha_grid))
    for (idx in seq_along(alpha_grid)) {
      alpha_val <- alpha_grid[idx]
      training_result <- train_model_cv(
        data_split$train,
        alpha_val,
        precomputed_results = train_precomputed,
        penalty_rate = 1.0,
        em_max_iter = em_max_iter,
        em_eps = em_eps
      )
      
      test_kl_loss <- evaluate_model_cv(training_result, data_split$test, alpha_val)
      num_edges <- sum(training_result$model$adjacency_matrix == 1)
      
      fold_rows[[idx]] <- data.frame(
        alpha = alpha_val,
        fold = fold,
        test_kl_loss = test_kl_loss,
        num_edges = num_edges,
        stringsAsFactors = FALSE
      )
    }
    
    do.call(rbind, fold_rows)
  }
  
  if (use_parallel && n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      # PSOCK workers (Windows / no fork): load the installed package so the
      # compiled C++ and internal R functions are available in each worker,
      # then ship the fold worker (its environment carries the data).
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, {
        suppressMessages(requireNamespace("ZICDT", quietly = TRUE))
      })
      parallel::clusterExport(cl, varlist = "process_single_fold",
                              envir = environment())
      fold_results <- parallel::parLapply(cl, folds, process_single_fold)
    } else {
      # fork-based workers (Unix/macOS) inherit the loaded namespace directly.
      fold_results <- parallel::mclapply(folds, process_single_fold,
                                         mc.cores = n_cores)
    }
  } else {
    fold_results <- lapply(folds, process_single_fold)
  }
  
  cv_results <- do.call(rbind, fold_results)
  rownames(cv_results) <- NULL
  cv_results <- cv_results[order(cv_results$fold, cv_results$alpha), , drop = FALSE]
  return(cv_results)
}

# Select optimal alpha based on CV results
# @param cv_results Results from cross_validate_alpha
# @param tie_breaking If TRUE, prefer simpler models (fewer edges) in case of ties
# @return List with optimal alpha and summary statistics
select_optimal_alpha <- function(cv_results, tie_breaking = TRUE) {
  # Compute summary statistics for each alpha using safer approach
  alpha_summary <- aggregate(test_kl_loss ~ alpha, cv_results, 
                             function(x) c(mean = mean(x), sd = sd(x)))
  alpha_summary$num_edges_mean <- aggregate(num_edges ~ alpha, cv_results, mean)$num_edges
  alpha_summary$num_edges_sd <- aggregate(num_edges ~ alpha, cv_results, sd)$num_edges
  
  # Extract mean scores (handling the matrix structure from aggregate)
  if (is.matrix(alpha_summary$test_kl_loss)) {
    mean_kl_losses <- alpha_summary$test_kl_loss[, "mean"]
  } else {
    mean_kl_losses <- alpha_summary$test_kl_loss
  }
  
  min_idx <- which.min(mean_kl_losses)
  min_kl_loss <- mean_kl_losses[min_idx]
  
  # Handle ties by preferring simpler models (fewer edges)
  if (tie_breaking) {
    # Find all alphas within tolerance of minimum KL loss
    tolerance <- 1e-6
    candidate_indices <- which(abs(mean_kl_losses - min_kl_loss) < tolerance)
    
    if (length(candidate_indices) > 1) {
      # Among candidates, choose the one with fewest edges
      candidate_edges <- alpha_summary$num_edges_mean[candidate_indices]
      optimal_idx <- candidate_indices[which.min(candidate_edges)]
    } else {
      optimal_idx <- min_idx
    }
  } else {
    optimal_idx <- min_idx
  }
  
  optimal_alpha <- alpha_summary$alpha[optimal_idx]
  
  cat("================================================================================\n")
  cat("ALPHA SELECTION RESULTS\n")
  cat("================================================================================\n")
  cat(sprintf("Optimal alpha: %.2f\n", optimal_alpha))
  
  # Extract values safely
  if (is.matrix(alpha_summary$test_kl_loss)) {
    mean_kl_loss <- alpha_summary$test_kl_loss[optimal_idx, "mean"]
    sd_kl_loss <- alpha_summary$test_kl_loss[optimal_idx, "sd"]
  } else {
    mean_kl_loss <- alpha_summary$test_kl_loss[optimal_idx]
    sd_kl_loss <- NA_real_
  }
  
  cat(sprintf("Mean CV KL loss: %.6f +/- %.6f\n", mean_kl_loss, sd_kl_loss))
  cat(sprintf("Mean edges: %.1f +/- %.1f\n",
              alpha_summary$num_edges_mean[optimal_idx],
              alpha_summary$num_edges_sd[optimal_idx]))
  cat("================================================================================\n")
  
  return(list(
    optimal_alpha = optimal_alpha,
    mean_kl_loss = mean_kl_loss,
    sd_kl_loss = sd_kl_loss,
    mean_edges = alpha_summary$num_edges_mean[optimal_idx],
    sd_edges = alpha_summary$num_edges_sd[optimal_idx],
    summary = alpha_summary,
    cv_results = cv_results
  ))
}

# Train final model with optimal alpha
# @param data_list Full dataset
# @param optimal_alpha Selected alpha value
# @param em_max_iter EM maximum iterations (default: 20000)
# @param em_eps EM convergence tolerance (default: 1e-8)
# @return List containing final model and learned parameters
train_final_model <- function(data_list,
                              optimal_alpha,
                              em_max_iter = 20000,
                              em_eps = 1e-8) {
  cat("================================================================================\n")
  cat("TRAINING FINAL MODEL WITH OPTIMAL ALPHA\n")
  cat("================================================================================\n")
  cat(sprintf("Training on full dataset with alpha = %.2f\n", optimal_alpha))
  
  # Train on FULL dataset with optimal alpha
  final_model <- edmonds_forest(data_list,
                                penalty_rate = 1.0,    # Fixed at 1.0
                                alpha = optimal_alpha,  # Use optimal alpha
                                em_max_iter = em_max_iter,
                                em_eps = em_eps)
  
  # Extract learned parameters for each node
  G_learned <- final_model$adjacency_matrix
  learned_params <- list()
  
  for (node in seq_len(nrow(G_learned))) {
    parent_ids <- which(G_learned[, node] == 1)
    
    if (length(parent_ids) == 0) {
      # No parents: use empirical mean
      y <- as.matrix(data_list[[node]])
      eta <- colMeans(y)
      learned_params[[node]] <- list(eta = eta, type = "no_parents")
    } else {
      # Has parents: use EM parameters
      parent_data <- data_list[parent_ids]
      y <- as.matrix(data_list[[node]])
      em_result <- EMalgorithm_cpp(parent_data, y, em_max_iter, em_eps)
      learned_params[[node]] <- list(
        w = em_result$w,
        eta = em_result$eta, 
        M = em_result$M,
        type = "with_parents"
      )
    }
  }
  
  cat(sprintf("Final model: %d edges, score = %.6f\n", 
              sum(final_model$adjacency_matrix == 1), final_model$total_score))
  cat("================================================================================\n")
  
  # Return both model and learned parameters
  return(list(
    model = final_model,
    learned_params = learned_params
  ))
}

# Complete cross-validation workflow
# @param data_list List of p matrices, each N x q_j
# @param alpha_grid Vector of alpha values to test (default: 0.01 to 1.0)
# @param K Number of CV folds (default: 5)
# @param tie_breaking If TRUE, prefer simpler models in case of ties
# @param use_parallel Whether to use parallel processing (default: TRUE)
# @param n_cores Number of cores to use (default: detectCores() - 1)
# @return Complete CV results and final model
cross_validate_alpha_complete <- function(data_list, 
                                          alpha_grid = NULL, 
                                          K = 10, 
                                          tie_breaking = TRUE,
                                          use_parallel = TRUE,
                                          n_cores = NULL,
                                          precompute_parallel = FALSE,
                                          precompute_n_cores = NULL,
                                          em_max_iter = 20000,
                                          em_eps = 1e-8) {
  # Default alpha grid if not provided
  if (is.null(alpha_grid)) {
    # Dense grid on (0, 2]: 2000 points, step 0.001. Reaching an interior optimum
    # instead of clamping is cheap here because each alpha reuses the cached EM fits.
    alpha_grid <- seq(0.001, 2, by = 0.001)
  }
  
  # Step 1: Cross-validation
  cv_results_df <- cross_validate_alpha(
    data_list = data_list,
    alpha_grid = alpha_grid,
    K = K,
    use_parallel = use_parallel,
    n_cores = n_cores,
    precompute_parallel = precompute_parallel,
    precompute_n_cores = precompute_n_cores,
    em_max_iter = em_max_iter,
    em_eps = em_eps
  )
  
  # Step 2: Select optimal alpha
  alpha_selection <- select_optimal_alpha(cv_results_df, tie_breaking)
  
  # Step 3: Train final model
  final_result <- train_final_model(
    data_list,
    alpha_selection$optimal_alpha,
    em_max_iter = em_max_iter,
    em_eps = em_eps
  )
  final_model <- final_result$model
  final_learned_params <- final_result$learned_params
  
  # Compile complete results
  complete_results <- list(
    # Cross-validation results
    cv_results = cv_results_df,
    cv_summary = cv_results_df,
    alpha_selection = alpha_selection,
    
    # Final model
    final_model = final_model,
    final_structure = final_model$adjacency_matrix,
    final_learned_params = final_learned_params,  # NEW: Learned parameters
    final_score = final_model$total_score,
    final_train_kl_loss = (final_model$total_score - 1.0 * sum(final_model$adjacency_matrix == 1)) / alpha_selection$optimal_alpha,
    final_num_edges = sum(final_model$adjacency_matrix == 1),
    
    # Model details
    node_names = rownames(final_model$adjacency_matrix),
    n_edges = final_model$n_edges,
    n_roots = final_model$n_roots,
    root_names = final_model$root_names,
    is_tree = final_model$is_tree,
    is_forest = final_model$is_forest,
    runtime = final_model$runtime,
    edmonds_cost = final_model$edmonds_cost,
    score_diff = final_model$score_diff
  )
  
  return(complete_results)
}

# ================================================================================
# SAVE RESULTS TO TEXT FILE WITH DESCRIPTIVE NAMING
# ================================================================================
# Save Edmonds algorithm results to a text file with descriptive naming
# @param forest_id String identifier for the forest (e.g., "forest_121_0000000100010000")
# @param dataset_id String identifier for the dataset (e.g., "dat_1")
# @param data_list List of node data matrices
# @param alpha_grid Vector of alpha values to test
# @param K Number of CV folds (default: 5)
# @param tie_breaking Whether to prefer simpler models in case of ties (default: TRUE)
# @param use_parallel Whether to use parallel processing (default: TRUE)
# @param n_cores Number of cores to use (default: detectCores() - 1)
# @param output_dir Directory to save results (default: current directory)
# @param precompute_parallel Whether to parallelize node precomputation within each fold
# @param precompute_n_cores Number of cores for node-level precomputation (NULL to auto-select)
# @param em_max_iter EM maximum iterations for both CV and final training
# @param em_eps EM convergence tolerance for both CV and final training
# @return List containing text file path, RData file path, and complete R results
