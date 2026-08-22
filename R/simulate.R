# Clean data simulator for ZICDT (replaces the driver-tangled generator in the
# original research script). Matches the paper's generative model.

#' Draw rows from a Dirichlet distribution
#'
#' @param n Number of rows to draw.
#' @param alpha Concentration vector (length \code{d}).
#' @return An \code{n x d} matrix; each row sums to 1.
#' @keywords internal
#' @noRd
rdirichlet_rows <- function(n, alpha) {
  d <- length(alpha)
  g <- matrix(stats::rgamma(n * d, shape = alpha, rate = 1), nrow = n, byrow = TRUE)
  g / rowSums(g)
}

#' Simulate zero-inflated compositional data on a directed tree
#'
#' Generates data from the paper's model: for each edge \eqn{k \to j} the child
#' conditional mean is \eqn{\omega_0 \eta_j + \omega_1 M_{jk} x^{(k)}} with a
#' column-stochastic transition matrix \eqn{M_{jk}}; observations are drawn as
#' \eqn{x^{(j)} \sim \mathrm{Dir}(\kappa \cdot \widehat{x}^{(j)})}.
#'
#' @param true_G A \code{p x p} 0/1 adjacency matrix of the true directed tree,
#'   where \code{true_G[k, j] == 1} means \eqn{k \to j} (k is the parent of j).
#'   Each column must have at most one non-zero entry (at most one parent).
#' @param d Integer vector of length \code{p} giving each node's compositional
#'   dimension, or a single integer used for all nodes.
#' @param n Number of samples (subjects).
#' @param omega1 Parental-influence weight (\eqn{\omega_1}); \eqn{\omega_0 = 1 - \omega_1}.
#' @param precision Dirichlet precision \eqn{\kappa} controlling dispersion.
#' @param zero_inflate Optional fraction in [0, 1); if > 0, that fraction of
#'   entries in each row is randomly set to exactly zero and the row renormalized,
#'   to mimic zero inflation.
#' @param seed Optional integer seed.
#' @return A list with \code{data} (a length-\code{p} list of \code{n x d_j}
#'   matrices) and \code{true_G} (the adjacency matrix used).
#' @export
zicdt_simulate <- function(true_G, d = 50, n = 100, omega1 = 0.95,
                           precision = 20, zero_inflate = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- nrow(true_G)
  if (length(d) == 1L) d <- rep(as.integer(d), p)
  if (length(d) != p) stop("`d` must have length 1 or nrow(true_G).")
  if (any(colSums(true_G != 0) > 1)) stop("Each node may have at most one parent (>=1 column of true_G has >1 entry).")

  omega0 <- 1 - omega1
  # topological order so parents are simulated before children
  ord <- topo_order(true_G)
  # baseline eta: root nodes are drawn flat, Dir(1); non-root nodes Dir(1/d).
  # (A spiky root would otherwise propagate through the tree and over-inflate zeros.)
  is_root <- colSums(true_G != 0) == 0
  eta <- lapply(seq_len(p), function(j) {
    conc <- if (is_root[j]) rep(1, d[j]) else rep(1 / d[j], d[j])
    as.numeric(rdirichlet_rows(1, conc))
  })
  Mlist <- vector("list", p)
  data <- vector("list", p)

  for (j in ord) {
    parent <- which(true_G[, j] != 0)
    if (length(parent) == 0) {
      xhat <- matrix(eta[[j]], nrow = n, ncol = d[j], byrow = TRUE)
    } else {
      k <- parent[1]
      # column-stochastic M: each of the d[k] columns is a composition on the child simplex
      M <- vapply(seq_len(d[k]), function(.) as.numeric(rdirichlet_rows(1, rep(1 / d[j], d[j]))),
                  numeric(d[j]))               # d[j] x d[k]
      Mlist[[j]] <- M
      parent_driven <- data[[k]] %*% t(M)      # n x d[j]
      xhat <- omega0 * matrix(eta[[j]], nrow = n, ncol = d[j], byrow = TRUE) +
        omega1 * parent_driven
    }
    xhat <- xhat / rowSums(xhat)
    X <- t(apply(xhat, 1, function(m) as.numeric(rdirichlet_rows(1, precision * m))))
    if (zero_inflate > 0) X <- inject_zeros(X, zero_inflate)
    data[[j]] <- X
  }
  names(data) <- as.character(seq_len(p))
  list(data = data, true_G = true_G, eta = eta, M = Mlist)
}

# Topological order of a directed tree/forest (parents before children).
topo_order <- function(G) {
  p <- nrow(G)
  indeg <- colSums(G != 0)
  ord <- integer(0); ready <- which(indeg == 0)
  while (length(ready) > 0) {
    v <- ready[1]; ready <- ready[-1]; ord <- c(ord, v)
    kids <- which(G[v, ] != 0)
    for (w in kids) { indeg[w] <- indeg[w] - 1; if (indeg[w] == 0) ready <- c(ready, w) }
  }
  if (length(ord) != p) stop("true_G is not acyclic.")
  ord
}

# Randomly set a fraction of entries per row to exactly zero, then renormalize.
inject_zeros <- function(X, rate) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    d <- ncol(X)
    k <- floor(rate * d)
    if (k > 0) {
      idx <- sample.int(d, k)
      X[i, idx] <- 0
      s <- sum(X[i, ])
      if (s > 0) X[i, ] <- X[i, ] / s
    }
  }
  X
}
