# Select the penalty by cross-validation and fit the final structure

Runs K-fold cross-validation over a grid of penalties, picks the penalty
that minimizes out-of-sample KL loss (with tie-breaking toward sparser
trees), and refits the final directed tree/forest on all data at the
selected penalty.

## Usage

``` r
zicdt_cv(
  data,
  alpha_grid = NULL,
  K = 5,
  em_max_iter = 20000,
  em_eps = 1e-08,
  use_parallel = FALSE,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  A named list of `n x d_j` compositional matrices (see
  [`zicdt_fit`](https://shuang-jie.github.io/ZICDT/reference/zicdt_fit.md)).

- alpha_grid:

  Numeric vector of candidate penalties. If `NULL`, a default grid
  `c(seq(0.01,0.09,0.01), seq(0.1,1,0.1))` is used.

- K:

  Number of cross-validation folds (default 5).

- em_max_iter, em_eps:

  EM iteration cap and convergence tolerance.

- use_parallel:

  Parallelize folds (default `FALSE`; sequential is the portable,
  package-safe path).

- seed:

  Optional integer. If supplied, the cross-validation fold split (the
  only stochastic step; edge EM fitting is deterministic) is seeded so
  the result is exactly reproducible. Pass a per-replicate seed to make
  a whole simulation study reproducible, including under parallel
  execution.

- verbose:

  If `FALSE` (default), suppress the engine's console output.

## Value

An object of class `"zicdt"` as in
[`zicdt_fit`](https://shuang-jie.github.io/ZICDT/reference/zicdt_fit.md),
with additional elements `cv` (the per-penalty CV table) and `params`
(the learned edge parameters of the final model).
