context("MVNormal C++ MCMC Implementation Tests")

test_that("MVNormal C++ MCMC produces valid results", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Generate test data
  set.seed(123)
  n <- 100
  d <- 2

  # True parameters for two clusters
  mu1 <- c(-2, -2)
  mu2 <- c(2, 2)
  Sigma1 <- diag(2) * 0.5
  Sigma2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)

  # Generate data
  data1 <- mvtnorm::rmvnorm(50, mu1, Sigma1)
  data2 <- mvtnorm::rmvnorm(50, mu2, Sigma2)
  y <- rbind(data1, data2)

  # Set up prior parameters
  g0Priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2),
    kappa0 = 1,
    nu = d + 2
  )

  # Create DP object
  dp <- DirichletProcessMvnormal(y, g0Priors)

  # Test with C++ backend
  set_use_cpp(TRUE)
  skip_if(!can_use_cpp(dp), "C++ backend not available for MVNormal")

  # Fit model
  dp_cpp <- Fit(dp, 100, progressBar = FALSE)

  # Basic validity checks
  expect_true(dp_cpp$numberClusters >= 1)
  expect_true(dp_cpp$numberClusters <= n)
  expect_equal(length(dp_cpp$clusterLabels), n)
  expect_true(all(dp_cpp$clusterLabels > 0))
  expect_equal(sum(dp_cpp$pointsPerCluster), n)

  # Check parameter dimensions
  expect_equal(dim(dp_cpp$clusterParameters$mu)[1], 1)
  expect_equal(dim(dp_cpp$clusterParameters$mu)[2], d)
  expect_equal(dim(dp_cpp$clusterParameters$sig)[1], d)
  expect_equal(dim(dp_cpp$clusterParameters$sig)[2], d)

  # Check chains exist
  expect_equal(length(dp_cpp$alphaChain), 100)
  expect_equal(length(dp_cpp$clusterLabelChain), 100)
})

test_that("MVNormal C++ vs R implementation equivalence", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Small dataset for comparison
  set.seed(456)
  y <- mvtnorm::rmvnorm(30, c(0, 0), diag(2))

  g0Priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2) * 2,
    kappa0 = 0.5,
    nu = 4
  )

  # Test R backend
  set_use_cpp(FALSE)
  set.seed(789)
  dp_r <- DirichletProcessMvnormal(y, g0Priors)
  dp_r <- Fit(dp_r, 50, progressBar = FALSE)

  # Test C++ backend
  set_use_cpp(TRUE)
  skip_if(!can_use_cpp(DirichletProcessMvnormal(y, g0Priors)))

  set.seed(789)
  dp_cpp <- DirichletProcessMvnormal(y, g0Priors)
  dp_cpp <- Fit(dp_cpp, 50, progressBar = FALSE)

  # Compare number of clusters (allowing for MCMC variability)
  expect_true(abs(dp_r$numberClusters - dp_cpp$numberClusters) <= 3)

  # Both should have valid alpha values
  expect_true(dp_r$alpha > 0)
  expect_true(dp_cpp$alpha > 0)
})

test_that("MVNormal handles edge cases correctly", {
  # Test with 1D data (degenerate case)
  y_1d <- matrix(rnorm(20), ncol = 1)
  g0_1d <- list(
    mu0 = 0,
    Lambda = matrix(1, 1, 1),
    kappa0 = 1,
    nu = 2
  )

  set_use_cpp(TRUE)
  dp_1d <- DirichletProcessMvnormal(y_1d, g0_1d)
  skip_if(!can_use_cpp(dp_1d))

  expect_silent(dp_1d <- Fit(dp_1d, 10, progressBar = FALSE))
  expect_true(dp_1d$numberClusters >= 1)

  # Test with perfect clustering
  y_perfect <- rbind(
    matrix(rep(c(-5, -5), 10), ncol = 2, byrow = TRUE),
    matrix(rep(c(5, 5), 10), ncol = 2, byrow = TRUE)
  )

  dp_perfect <- DirichletProcessMvnormal(y_perfect)
  dp_perfect <- Fit(dp_perfect, 50, progressBar = FALSE)

  # Should find approximately 2 clusters
  expect_true(dp_perfect$numberClusters >= 1)
  expect_true(dp_perfect$numberClusters <= 5)
})
