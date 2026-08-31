# Competitor methods for the simulation study: PC-stable and LiNGAM, adapted
# (block-level) to compositional nodes. Each function takes `blocks` (a list of
# n x d_j matrices) and returns a p x p directed forest adjacency:
# forest[i, j] == 1 means edge i -> j.
#
# PC-stable is provided in two variants:
#   pc_block_forest      : PC on the raw (scaled) coordinates.
#   pc_ilr_block_forest  : PC on the ILR (isometric log-ratio) coordinates, the
#                          statistically correct representation for compositions.
# Both default to conservative = FALSE and u2pd = "retry":
#   - conservative = FALSE lets PC orient v-structures + Meek rules (conservative
#     = TRUE leaves most edges unoriented, which the block vote then breaks
#     arbitrarily -> poor, order-dependent direction recovery).
#   - u2pd = "retry" returns a coherent, extendable DAG, avoiding a non-extendable
#     PDAG whose downstream repair would inject index-order information.
# ILR uses logs, so it is undefined on exact zeros; pc_ilr_block_forest therefore
# replaces zeros first (multiplicative replacement by default). This imputation is
# an intrinsic requirement of log-ratio methods on zero-inflated data -- the model
# in this package needs no such step.

suppressMessages({library(pcalg)})

## ---- shared block-level majority vote: coordinate CPDAG -> directed forest ----
.block_counts <- function(A, bs) {
  nb <- length(bs); bstart <- cumsum(c(1, bs[-nb])); bend <- cumsum(bs)
  counts <- matrix(0, nb, nb)
  for (i in seq_len(nb)) for (j in seq_len(nb)) {
    if (i == j) next
    counts[i, j] <- sum(A[bstart[i]:bend[i], bstart[j]:bend[j]] == 1)
  }
  counts
}
.forest_from_counts <- function(counts) {
  nb <- nrow(counts); dagm <- matrix(0, nb, nb)
  # pick the majority direction per block pair; ties (rare with u2pd="retry")
  # fall back to i<j purely to avoid dropping a detected edge.
  for (i in seq_len(nb - 1)) for (j in (i + 1):nb) {
    if (counts[i, j] > counts[j, i]) dagm[i, j] <- 1
    else if (counts[j, i] > counts[i, j]) dagm[j, i] <- 1
    else if (counts[i, j] > 0) dagm[i, j] <- 1
  }
  forest <- matrix(0, nb, nb)
  for (child in seq_len(nb)) {
    parents <- which(dagm[, child] == 1)
    if (length(parents) == 0) next
    forest[parents[which.max(counts[parents, child])], child] <- 1
  }
  forest
}

## ---- PC-stable on raw (scaled) coordinates ----
pc_block_forest <- function(blocks, alpha = 0.05, conservative = FALSE, u2pd = "retry") {
  U <- do.call(cbind, lapply(blocks, scale))
  colnames(U) <- as.character(seq_len(ncol(U)))
  pc_fit <- pc(suffStat = list(C = cor(U), n = nrow(U)), indepTest = gaussCItest,
               alpha = alpha, labels = colnames(U), skel.method = "stable",
               conservative = conservative, solve.confl = FALSE, u2pd = u2pd)
  A <- as(pc_fit@graph, "matrix")
  .forest_from_counts(.block_counts(A, sapply(blocks, ncol)))
}

## ---- zero replacement (needed before ILR: log is undefined at 0) ----
# Multiplicative simple replacement: zeros -> 0.65 * (smallest positive in row),
# non-zeros scaled down so the row still sums to 1.
.zero_replace_mult <- function(X) {
  X <- as.matrix(X)
  for (i in seq_len(nrow(X))) {
    xi <- X[i, ]; z <- xi == 0
    if (any(z) && any(!z)) {
      delta <- min(xi[xi > 0]) * 0.65
      xi[z] <- delta; xi[!z] <- xi[!z] * (1 - sum(z) * delta)
      xi <- xi / sum(xi)
    }
    X[i, ] <- xi
  }
  X
}

## ---- PC-stable on ILR coordinates (compositional-correct) ----
pc_ilr_block_forest <- function(blocks, alpha = 0.05, conservative = FALSE,
                                u2pd = "retry", zero_replace = .zero_replace_mult) {
  if (!requireNamespace("compositions", quietly = TRUE))
    stop("pc_ilr_block_forest needs the 'compositions' package.")
  reps <- lapply(blocks, function(X)
    as.matrix(compositions::ilr(compositions::acomp(zero_replace(X)))))
  bs <- sapply(reps, ncol)                 # each block -> d_j - 1 ILR coords
  U <- do.call(cbind, reps); colnames(U) <- as.character(seq_len(ncol(U)))
  pc_fit <- pc(suffStat = list(C = cor(U), n = nrow(U)), indepTest = gaussCItest,
               alpha = alpha, labels = colnames(U), skel.method = "stable",
               conservative = conservative, solve.confl = FALSE, u2pd = u2pd)
  A <- as(pc_fit@graph, "matrix")
  .forest_from_counts(.block_counts(A, bs))
}

## ---- LiNGAM (per-block PCA representation, then DirectLiNGAM) ----
.make_block_rep <- function(blocks, pca_var = 0.90, pca_max_per_block = NULL) {
  n <- nrow(blocks[[1]]); q <- length(blocks)
  if (is.null(pca_max_per_block)) pca_max_per_block <- max(1L, min(floor(n / q), 20L) - 1L)
  reps <- vector("list", q); col_names <- character(0)
  for (b in seq_len(q)) {
    Xb <- scale(blocks[[b]], center = TRUE, scale = FALSE)
    pc <- prcomp(Xb)
    var_prop <- summary(pc)$importance[2, ]; cumvar <- cumsum(var_prop)
    if (is.null(pca_var) || pca_var <= 0) k <- 1L else {
      k <- which(cumvar >= pca_var)[1]; if (is.na(k)) k <- length(cumvar)
    }
    k <- min(k, pca_max_per_block)
    scores <- pc$x[, seq_len(k), drop = FALSE]
    if (cor(scores[, 1], rowMeans(Xb)) < 0) scores[, 1] <- -scores[, 1]
    reps[[b]] <- scores
    col_names <- c(col_names, paste0("Block_", b, "_PC", seq_len(k)))
  }
  Z <- do.call(cbind, reps); colnames(Z) <- col_names; Z
}

lingam_block_forest <- function(blocks, pca_var = 0.90) {
  Z <- .make_block_rep(blocks, pca_var = pca_var)
  Z <- scale(Z, center = TRUE, scale = TRUE)
  fit <- lingam(Z)
  B <- fit$Bpruned; B[abs(B) < 1e-8] <- 0            # B[i, j] != 0 means j -> i
  rownames(B) <- colnames(B) <- colnames(Z)
  block_of <- as.integer(sub("Block_(\\d+)_PC.*", "\\1", colnames(Z)))
  nb <- max(block_of)
  counts <- matrix(0, nb, nb); p <- nrow(B)
  for (child in seq_len(p)) {
    parents <- which(B[child, ] != 0)
    if (length(parents) == 0) next
    best <- parents[which.max(abs(B[child, parents]))]  # PC-level best parent
    i <- block_of[best]; j <- block_of[child]
    if (i != j) counts[i, j] <- counts[i, j] + 1
  }
  .forest_from_counts(counts)
}

## ---- directed-edge recovery metric (TPR, FDR, MCC), shared by all methods ----
edge_metrics <- function(true_G, learned_G) {
  p <- nrow(true_G); TP <- FP <- FN <- TN <- 0
  for (i in 1:p) for (j in 1:p) {
    if (i == j) next
    t <- true_G[i, j] != 0; l <- learned_G[i, j] != 0
    if (t && l) TP <- TP + 1 else if (!t && l) FP <- FP + 1
    else if (t && !l) FN <- FN + 1 else TN <- TN + 1
  }
  TPR <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  FDR <- if (TP + FP > 0) FP / (TP + FP) else 0
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- if (den > 0) (TP*TN - FP*FN)/den else 0
  c(TPR = TPR, FDR = FDR, MCC = MCC)
}
