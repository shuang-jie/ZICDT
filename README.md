# ZICDT

**Structure Learning for Directed Trees with Zero-Inflated Compositional Nodes**

`ZICDT` learns a directed tree (or forest) in which **each node is an entire
zero-inflated composition** on the probability simplex. It scores candidate edges
with a Kullback–Leibler divergence (finite on exact zeros), models each child's
conditional mean as a mixture of a baseline and a parent-driven, column-stochastic
transition map, fits edge parameters by an EM algorithm (C++ via Rcpp), finds the
globally optimal spanning arborescence with the **Chu–Liu/Edmonds** algorithm
(a virtual root accommodates forests), and selects the edge penalty by
cross-validation.

## Installation

```r
# from the package directory
devtools::install("ZICDT")     # or: R CMD INSTALL ZICDT
```

Requires a C++ compiler with `Rcpp` and `RcppArmadillo`.

## Quick start

```r
library(ZICDT)

## simulate a 5-node tree:  1->2->3 and 1->4->5
G <- matrix(0, 5, 5); G[1,2] <- 1; G[2,3] <- 1; G[1,4] <- 1; G[4,5] <- 1
sim <- zicdt_simulate(G, d = c(8,10,12,9,11), n = 100, seed = 7)

## fit at a fixed penalty
fit <- zicdt_fit(sim$data, alpha = 2)
print(fit)
zicdt_edges(fit)

## select the penalty by cross-validation
cvfit <- zicdt_cv(sim$data, alpha_grid = c(seq(0.1, 1, 0.1), 2, 5), K = 5)
print(cvfit)
cvfit$alpha    # selected penalty
cvfit$cv       # per-penalty CV table
```

## Main functions

| function | purpose |
|---|---|
| `zicdt_fit(data, alpha)` | learn the tree/forest at a fixed penalty |
| `zicdt_cv(data, alpha_grid, K)` | choose the penalty by K-fold CV and refit |
| `zicdt_simulate(true_G, d, n)` | simulate compositional data on a known tree |
| `zicdt_edges(object)` | learned edges as a `parent -> child` data frame |

`data` is a **named list** of `n x d_j` matrices — one per node — each row a
composition (nonnegative, summing to 1). Larger `alpha` weights the KL fit more
heavily and yields **denser** structures.

## Notes

- The C++ engine (`src/EM_DAG.cpp`) is the validated implementation from the paper.
- Cross-validation runs sequentially by default (`use_parallel = FALSE`) for
  portability; fork-based parallelism is available on Unix via `use_parallel = TRUE`.
