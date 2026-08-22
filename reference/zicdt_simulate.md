# Simulate zero-inflated compositional data on a directed tree

Generates data from the paper's model: for each edge \\k \to j\\ the
child conditional mean is \\\omega_0 \eta_j + \omega_1 M\_{jk} x^{(k)}\\
with a column-stochastic transition matrix \\M\_{jk}\\; observations are
drawn as \\x^{(j)} \sim \mathrm{Dir}(\kappa \cdot \widehat{x}^{(j)})\\.

## Usage

``` r
zicdt_simulate(
  true_G,
  d = 50,
  n = 100,
  omega1 = 0.95,
  precision = 20,
  zero_inflate = 0,
  seed = NULL
)
```

## Arguments

- true_G:

  A `p x p` 0/1 adjacency matrix of the true directed tree, where
  `true_G[k, j] == 1` means \\k \to j\\ (k is the parent of j). Each
  column must have at most one non-zero entry (at most one parent).

- d:

  Integer vector of length `p` giving each node's compositional
  dimension, or a single integer used for all nodes.

- n:

  Number of samples (subjects).

- omega1:

  Parental-influence weight (\\\omega_1\\); \\\omega_0 = 1 - \omega_1\\.

- precision:

  Dirichlet precision \\\kappa\\ controlling dispersion.

- zero_inflate:

  Optional fraction in \[0, 1); if \> 0, that fraction of entries in
  each row is randomly set to exactly zero and the row renormalized, to
  mimic zero inflation.

- seed:

  Optional integer seed.

## Value

A list with `data` (a length-`p` list of `n x d_j` matrices) and
`true_G` (the adjacency matrix used).
