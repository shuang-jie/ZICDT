# Getting started with ZICDT

## Overview

`ZICDT` learns a **directed tree (or forest)** in which each node is an
entire **zero-inflated composition** on the probability simplex. It
scores candidate edges with a Kullback–Leibler divergence (finite on
exact zeros), models each child’s conditional mean as a mixture of a
baseline and a parent-driven, column-stochastic transition map, fits
edge parameters by EM (in C++), searches for the globally optimal
spanning arborescence with the Chu–Liu/Edmonds algorithm, and selects
the edge penalty by cross-validation.

``` r

library(ZICDT)
```

## Simulate data on a known tree

Each node is an `n x d_j` matrix whose rows are compositions. Here we
build a 5-node tree with edges `1->2->3` and `1->4->5`.

``` r

G <- matrix(0, 5, 5)
G[1, 2] <- 1; G[2, 3] <- 1; G[1, 4] <- 1; G[4, 5] <- 1

sim <- zicdt_simulate(G, d = c(8, 10, 12, 9, 11), n = 100, seed = 7)
lapply(sim$data, dim)          # one n x d_j matrix per node
#> $`1`
#> [1] 100   8
#> 
#> $`2`
#> [1] 100  10
#> 
#> $`3`
#> [1] 100  12
#> 
#> $`4`
#> [1] 100   9
#> 
#> $`5`
#> [1] 100  11
```

## Fit at a fixed penalty

[`zicdt_fit()`](https://shuang-jie.github.io/ZICDT/reference/zicdt_fit.md)
learns the globally optimal structure for a given penalty `alpha`. Note
the parameterization: `alpha` multiplies the KL fit, so **larger `alpha`
yields a denser structure**.

``` r

fit <- zicdt_fit(sim$data, alpha = 2)
fit
#> ZICDT directed tree
#>   nodes: 5  |  edges: 4  |  roots: 1
#>   penalty (alpha): 2
#>   score: 117.9108
#>   edges:
#>     1 -> 2
#>     2 -> 3
#>     3 -> 4
#>     4 -> 5
zicdt_edges(fit)
#>   parent child
#> 1      1     2
#> 2      2     3
#> 3      3     4
#> 4      4     5
```

## Choose the penalty by cross-validation (main workflow)

[`zicdt_cv()`](https://shuang-jie.github.io/ZICDT/reference/zicdt_cv.md)
runs K-fold cross-validation over a grid of penalties, selects the one
minimizing out-of-sample KL loss, and refits the final tree. This is the
usual entry point.

``` r

cvfit <- zicdt_cv(sim$data, alpha_grid = c(seq(0.1, 1, 0.1), 2, 5), K = 3)
cvfit
#> ZICDT directed forest
#>   nodes: 5  |  edges: 2  |  roots: 3
#>   penalty (alpha): 0.5
#>   score: 31.6120
#>   edges:
#>     2 -> 3
#>     4 -> 5
cvfit$alpha        # selected penalty
#> [1] 0.5
head(cvfit$cv)     # per-penalty CV table
#>   alpha fold test_kl_loss num_edges
#> 1   0.1    1     26.52684         0
#> 2   0.2    1     26.52684         0
#> 3   0.3    1     25.19383         1
#> 4   0.4    1     25.19383         1
#> 5   0.5    1     24.37119         2
#> 6   0.6    1     24.37119         2
```

The learned directed tree is available as an adjacency matrix
(`cvfit$adjacency`, with `adjacency[k, j] == 1` meaning edge `k -> j`),
as a tidy edge list via `zicdt_edges(cvfit)`, and as an `igraph` object
(`cvfit$graph`). The fitted edge transition matrices are in
`cvfit$params`.
