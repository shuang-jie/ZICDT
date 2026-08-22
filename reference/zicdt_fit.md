# Fit a directed tree/forest at a fixed penalty

Learns the globally optimal directed tree (or forest) over compositional
nodes at a given edge penalty `alpha`, using EM edge fitting and the
Chu-Liu/Edmonds algorithm.

## Usage

``` r
zicdt_fit(
  data,
  alpha,
  penalty_rate = 1,
  em_max_iter = 5000,
  em_eps = 1e-06,
  verbose = FALSE
)
```

## Arguments

- data:

  A named list of `n x d_j` matrices, one per node; each row is a
  composition (nonnegative, sums to 1). All matrices must share the same
  number of rows `n`.

- alpha:

  Edge penalty (larger `alpha` -\> sparser structure).

- penalty_rate:

  Per-edge complexity penalty (default 1).

- em_max_iter, em_eps:

  EM iteration cap and convergence tolerance.

- verbose:

  If `FALSE` (default), suppress the engine's console output.

## Value

An object of class `"zicdt"`: a list with `adjacency` (\\k \to j\\ when
`adjacency[k, j] == 1`), `n_edges`, `n_roots`, `roots`, `is_tree`,
`score`, `alpha`, and `graph` (an igraph object).
