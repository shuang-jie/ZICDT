# Changelog

## ZICDT 0.1.8

- No-parent (root) loss: the fitted mean `eta = colMeans(y)` is now
  floored at 1e-100 in `precompute_node_results()` and in the C++
  `compute_node_marginal_cost`. A taxon whose only nonzero entries are
  subnormal (~1e-323, from Dirichlet underflow) had a mean of exactly 0,
  which made the no-parent loss infinite: the node could never be a
  root, and an edge into it was forced at every penalty. Results are
  unchanged whenever every taxon mean exceeds 1e-100.

## ZICDT 0.1.7

- EM: a starting `eta(r) <= 1e-100` (not only `eta(r) == 0`) is replaced
  by a small positive value and renormalised. A taxon with only
  subnormal nonzero entries otherwise gave 0/0 = NaN in the first
  E-step, and the EM never converged.

## ZICDT 0.1.6

- EM: a child taxon that is zero in every sample (`eta(r) = 0` at the
  start) no longer produces NaN; such entries start at a small positive
  value. Results are unchanged on data without such a taxon.
- New internal `EMalgorithm_init_cpp()`: the same EM from a
  caller-supplied start (`w`, `eta`, `M`), for multi-start diagnostics.

## ZICDT 0.1.2

- zicdt_cv() gains a `seed` argument that seeds the cross-validation
  fold split (the only stochastic step; EM fitting is deterministic),
  making results exactly reproducible, including under parallel
  execution.

## ZICDT 0.1.1

- zicdt_simulate(): draw root-node baselines from Dir(1) (non-root from
  Dir(1/d)), matching the paper generator. Fixes over-inflated zeros
  (~49%) that had collapsed structure recovery on simulated data.

## ZICDT 0.1.0

- Initial release.
- [`zicdt_cv()`](https://shuang-jie.github.io/ZICDT/reference/zicdt_cv.md)
  learns a directed tree/forest with the edge penalty chosen by K-fold
  cross-validation;
  [`zicdt_fit()`](https://shuang-jie.github.io/ZICDT/reference/zicdt_fit.md)
  fits at a fixed penalty.
- [`zicdt_simulate()`](https://shuang-jie.github.io/ZICDT/reference/zicdt_simulate.md)
  generates zero-inflated compositional data on a known tree.
- C++ EM engine via Rcpp/RcppArmadillo; Chu-Liu/Edmonds tree search.
