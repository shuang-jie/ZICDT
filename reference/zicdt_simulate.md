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
  no_dead_cols = FALSE,
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

  Optional target zero fraction in \[0, 1). If \> 0, each observation is
  forced to have exactly `floor(zero_inflate * d_j)` zero components:
  the smallest-abundance entries (including any that were already zero
  from Dirichlet underflow) become structural zeros and the row is
  renormalized. This yields a zero-inflated Dirichlet whose realized
  zero rate equals `zero_inflate` (provided it exceeds the small natural
  rate), with zeros occurring at low abundance as in real compositional
  data.

- no_dead_cols:

  Logical; if `TRUE`, after zero-injection any column (component) that
  became zero in every observation is rescued by restoring its largest
  original value in one subject, guaranteeing every component is present
  in at least one observation (as in real data). Default `FALSE`.

- seed:

  Optional integer seed.

## Value

A list with `data` (a length-`p` list of `n x d_j` matrices) and
`true_G` (the adjacency matrix used).
