# ZICDT 0.1.1

* zicdt_simulate(): draw root-node baselines from Dir(1) (non-root from Dir(1/d)),
  matching the paper generator. Fixes over-inflated zeros (~49%) that had collapsed
  structure recovery on simulated data.

# ZICDT 0.1.0

* Initial release.
* `zicdt_cv()` learns a directed tree/forest with the edge penalty chosen by
  K-fold cross-validation; `zicdt_fit()` fits at a fixed penalty.
* `zicdt_simulate()` generates zero-inflated compositional data on a known tree.
* C++ EM engine via Rcpp/RcppArmadillo; Chu-Liu/Edmonds tree search.
