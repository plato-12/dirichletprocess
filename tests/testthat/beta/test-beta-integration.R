context("Beta DP Integration Tests")

test_that("Full Beta DP MCMC recovers known parameters", {
  set.seed(1234)
  # Generate data from known mixture
  n1 <- 50
  n2 <- 50
  data <- c(
    rbeta(n1, 2, 8),   # Low mean cluster
    rbeta(n2, 8, 2)    # High mean cluster
  )

  # Fit DP model
  dp <- DirichletProcessBeta(data, verbose = FALSE)
  dp <- Fit(dp, its = 500, progressBar = FALSE)

  # Check clustering
  expect_true(dp$numberClusters >= 2)
  expect_true(dp$numberClusters <= 4)  # Shouldn't overfit too much

  # Check parameter recovery
  mu_values <- sort(dp$clusterParameters[[1]])

  # Should have one cluster with low mean and one with high mean
  expect_true(any(mu_values < 0.3))  # Low cluster
  expect_true(any(mu_values > 0.7))  # High cluster
})

test_that("Beta DP handles single cluster data", {
  set.seed(2345)
  # Generate data from single Beta
  data <- rbeta(100, 5, 5)  # Symmetric around 0.5

  dp <- DirichletProcessBeta(data, verbose = FALSE)
  dp <- Fit(dp, its = 200, progressBar = FALSE)

  # Should mostly find 1 cluster
  expect_true(dp$numberClusters <= 3)

  # Parameter should be near true value
  if (dp$numberClusters == 1) {
    expect_equal(dp$clusterParameters[[1]][1], 0.5, tolerance = 0.1)
  }
})

test_that("Beta DP alpha updates work correctly", {
  set.seed(3456)
  data <- rbeta(50, 3, 3)

  dp <- DirichletProcessBeta(
    data,
    alphaPriors = c(2, 0.5),  # Prior mean = 4
    verbose = FALSE
  )

  initial_alpha <- dp$alpha

  # Run MCMC with alpha updates
  dp <- Fit(dp, its = 100, progressBar = FALSE, updatePrior = TRUE)

  # Alpha should have changed
  expect_true(dp$alpha != initial_alpha)

  # Should be reasonable value
  expect_true(dp$alpha > 0.1)
  expect_true(dp$alpha < 20)
})
