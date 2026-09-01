# User-facing API for ZICDT. Thin wrappers over the validated engine
# (edmonds.R + cv.R + src/EM_DAG.cpp).

#' Fit a directed tree/forest at a fixed penalty
#'
#' Learns the globally optimal directed tree (or forest) over compositional nodes
#' at a given edge penalty \code{alpha}, using EM edge fitting and the
#' Chu-Liu/Edmonds algorithm.
#'
#' @param data A named list of \code{n x d_j} matrices, one per node; each row is
#'   a composition (nonnegative, sums to 1). All matrices must share the same
#'   number of rows \code{n}.
#' @param alpha Edge penalty (larger \code{alpha} -> sparser structure).
#' @param penalty_rate Per-edge complexity penalty (default 1).
#' @param em_max_iter,em_eps EM iteration cap and convergence tolerance.
#' @param verbose If \code{FALSE} (default), suppress the engine's console output.
#' @return An object of class \code{"zicdt"}: a list with \code{adjacency}
#'   (\eqn{k \to j} when \code{adjacency[k, j] == 1}), \code{n_edges},
#'   \code{n_roots}, \code{roots}, \code{is_tree}, \code{score}, \code{alpha},
#'   and \code{graph} (an \pkg{igraph} object).
#' @export
zicdt_fit <- function(data, alpha, penalty_rate = 1,
                      em_max_iter = 5000, em_eps = 1e-6, verbose = FALSE) {
  stopifnot(is.list(data), length(data) >= 2)
  run <- function() edmonds_forest(data, penalty_rate = penalty_rate, alpha = alpha,
                                   em_max_iter = em_max_iter, em_eps = em_eps)
  res <- if (verbose) run() else { utils::capture.output(res <- run()); res }
  out <- list(
    adjacency = res$adjacency_matrix,
    n_edges   = res$n_edges,
    n_roots   = res$n_roots,
    roots     = res$root_names,
    is_tree   = res$is_tree,
    score     = res$total_score,
    alpha     = alpha,
    graph     = res$graph
  )
  class(out) <- "zicdt"
  out
}

#' Select the penalty by cross-validation and fit the final structure
#'
#' Runs K-fold cross-validation over a grid of penalties, picks the penalty that
#' minimizes out-of-sample KL loss (with tie-breaking toward sparser trees), and
#' refits the final directed tree/forest on all data at the selected penalty.
#'
#' @param data A named list of \code{n x d_j} compositional matrices (see \code{\link{zicdt_fit}}).
#' @param alpha_grid Numeric vector of candidate penalties. If \code{NULL}, a
#'   default grid \code{c(seq(0.01,0.09,0.01), seq(0.1,1,0.1))} is used.
#' @param K Number of cross-validation folds (default 5).
#' @param em_max_iter,em_eps EM iteration cap and convergence tolerance.
#' @param use_parallel Parallelize folds (default \code{FALSE}; sequential is the
#'   portable, package-safe path).
#' @param seed Optional integer. If supplied, the cross-validation fold split
#'   (the only stochastic step; edge EM fitting is deterministic) is seeded so
#'   the result is exactly reproducible. Pass a per-replicate seed to make a
#'   whole simulation study reproducible, including under parallel execution.
#' @param verbose If \code{FALSE} (default), suppress the engine's console output.
#' @return An object of class \code{"zicdt"} as in \code{\link{zicdt_fit}}, with
#'   additional elements \code{cv} (the per-penalty CV table) and \code{params}
#'   (the learned edge parameters of the final model).
#' @export
zicdt_cv <- function(data, alpha_grid = NULL, K = 10,
                     em_max_iter = 20000, em_eps = 1e-8,
                     use_parallel = FALSE, seed = NULL, verbose = FALSE) {
  stopifnot(is.list(data), length(data) >= 2)
  if (!is.null(seed)) set.seed(seed)
  run <- function() cross_validate_alpha_complete(
    data_list = data, alpha_grid = alpha_grid, K = K,
    use_parallel = use_parallel, em_max_iter = em_max_iter, em_eps = em_eps)
  res <- if (verbose) run() else { utils::capture.output(res <- run()); res }
  # Candidate-edge EM convergence on the full data: precompute_node_results calls
  # the same EMalgorithm_cpp used by the scorer, so em_iters here == the scorer's.
  # (Adds one full-data candidate-fit pass; em_iters == em_max_iter means it hit the cap.)
  emc <- tryCatch({
    pnr <- precompute_node_results(data, em_max_iter = em_max_iter, em_eps = em_eps)
    rows <- lapply(seq_along(pnr), function(j) {
      nr <- pnr[[j]]; pk <- setdiff(names(nr), "0")
      if (length(pk) == 0) return(NULL)
      data.frame(child = j, parent = as.integer(pk),
                 em_iters = vapply(pk, function(k) as.integer(nr[[k]]$em_iters), integer(1)),
                 row.names = NULL)
    })
    do.call(rbind, rows)
  }, error = function(e) NULL)
  out <- list(
    adjacency = res$final_structure,
    n_edges   = res$final_num_edges,
    n_roots   = res$n_roots,
    roots     = res$node_names[res$final_model$roots],
    is_tree   = isTRUE(res$final_model$is_tree),
    score     = res$final_score,
    alpha     = res$alpha_selection$optimal_alpha,
    cv        = res$cv_results,
    cv_curve  = tryCatch(stats::aggregate(test_kl_loss ~ alpha, res$cv_results, mean),
                         error = function(e) NULL),   # mean held-out KL loss per alpha
    em_convergence = emc,                             # per (child,parent) actual EM iters
    em_iters_max   = if (!is.null(emc)) max(emc$em_iters, na.rm = TRUE) else NA_integer_,
    em_n_hitcap    = if (!is.null(emc)) sum(emc$em_iters >= em_max_iter, na.rm = TRUE) else NA_integer_,
    em_max_iter    = em_max_iter,
    params    = res$final_learned_params,
    graph     = res$final_model$graph
  )
  class(out) <- "zicdt"
  out
}

#' Extract the learned edges as a data frame
#'
#' @param object A \code{"zicdt"} object.
#' @return A data frame with columns \code{parent} and \code{child}
#'   (node names), one row per directed edge \eqn{parent \to child}.
#' @export
zicdt_edges <- function(object) {
  A <- object$adjacency
  idx <- which(A != 0, arr.ind = TRUE)
  nm <- rownames(A); if (is.null(nm)) nm <- as.character(seq_len(nrow(A)))
  if (nrow(idx) == 0) return(data.frame(parent = character(0), child = character(0)))
  data.frame(parent = nm[idx[, 1]], child = nm[idx[, 2]], stringsAsFactors = FALSE)
}

#' @export
print.zicdt <- function(x, ...) {
  cat(sprintf("ZICDT directed %s\n", if (isTRUE(x$is_tree)) "tree" else "forest"))
  cat(sprintf("  nodes: %d  |  edges: %d  |  roots: %d\n",
              nrow(x$adjacency), x$n_edges, x$n_roots))
  if (!is.null(x$alpha)) cat(sprintf("  penalty (alpha): %.4g\n", x$alpha))
  if (!is.null(x$score)) cat(sprintf("  score: %.4f\n", x$score))
  e <- zicdt_edges(x)
  if (nrow(e) > 0) {
    cat("  edges:\n")
    for (i in seq_len(nrow(e))) cat(sprintf("    %s -> %s\n", e$parent[i], e$child[i]))
  }
  invisible(x)
}
