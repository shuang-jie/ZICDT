test_that("simulate returns valid compositions", {
  G <- matrix(0, 4, 4); G[1,2] <- 1; G[2,3] <- 1; G[1,4] <- 1
  sim <- zicdt_simulate(G, d = 6, n = 30, seed = 1)
  expect_length(sim$data, 4)
  expect_equal(nrow(sim$data[[1]]), 30)
  expect_true(all(abs(rowSums(sim$data[[2]]) - 1) < 1e-8))
})

test_that("zicdt_fit runs and returns a valid structure", {
  G <- matrix(0, 4, 4); G[1,2] <- 1; G[2,3] <- 1; G[1,4] <- 1
  sim <- zicdt_simulate(G, d = 6, n = 40, seed = 2)
  fit <- zicdt_fit(sim$data, alpha = 2)
  expect_s3_class(fit, "zicdt")
  expect_true(max(colSums(fit$adjacency)) <= 1)   # each node <= 1 parent
  expect_true(fit$n_edges >= 0)
})
