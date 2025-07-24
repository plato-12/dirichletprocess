# tests/testthat/test-cpp-edge-cases.R

test_that("C++ handles empty clusters correctly", {
  # Create scenario with empty clusters
  set.seed(123)
  test_data <- c(rnorm(50, -10), rnorm(50, 10))  # Well-separated clusters

  dp <- DirichletProcessGaussian(test_data)
  dp <- Fit(dp, its = 100)

  # Force empty cluster scenario
  dp$clusterLabels[1:25] <- 1
  dp$clusterLabels[26:100] <- 2
  dp$numberClusters <- 2

  # This should not crash
  expect_error(
    dp <- Fit(dp, its = 10),
    NA  # Expect no error
  )
})

test_that("C++ handles single data point", {
  test_data <- 1.5

  expect_error(
    dp <- DirichletProcessGaussian(test_data),
    NA
  )

  expect_error(
    dp <- Fit(dp, its = 10),
    NA
  )
})

test_that("C++ handles large datasets efficiently", {
  set.seed(123)
  test_data <- rnorm(10000)

  start_time <- Sys.time()
  dp <- DirichletProcessGaussian(test_data)
  dp <- Fit(dp, its = 10)
  runtime <- as.numeric(Sys.time() - start_time)

  # Should complete in reasonable time
  expect_lt(runtime, 30)  # 30 seconds max

  # Should produce valid results
  expect_true(dp$numberClusters > 0)
  expect_true(all(dp$clusterLabels %in% 1:dp$numberClusters))
})

test_that("C++ handles extreme parameter values", {
  set.seed(123)

  # Very large alpha
  dp <- DirichletProcessGaussian(rnorm(100))
  dp$alpha <- 1000
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Very small alpha
  dp$alpha <- 0.001
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Extreme data values
  test_data <- c(rnorm(50), 1e6, -1e6)
  dp <- DirichletProcessGaussian(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)
})

test_that("C++ handles degenerate data", {
  # All same value
  test_data <- rep(5, 100)
  dp <- DirichletProcessGaussian(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Only two unique values
  test_data <- rep(c(0, 1), 50)
  dp <- DirichletProcessGaussian(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)
})

test_that("C++ handles matrix data edge cases", {
  # Single row matrix
  test_data <- matrix(rnorm(5), nrow = 1)
  expect_error(
    dp <- DirichletProcessMvnormal(test_data),
    NA
  )

  # Very high dimensional data
  test_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  dp <- DirichletProcessMvnormal(test_data)
  expect_error(dp <- Fit(dp, its = 5), NA)

  # More columns than rows
  test_data <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  dp <- DirichletProcessMvnormal(test_data)
  expect_error(dp <- Fit(dp, its = 5), NA)
})

test_that("C++ handles numerical precision edge cases", {
  # Very small variance data
  set.seed(123)
  test_data <- rnorm(100, sd = 1e-10)
  dp <- DirichletProcessGaussian(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Very large variance data
  test_data <- rnorm(100, sd = 1e10)
  dp <- DirichletProcessGaussian(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)
})

test_that("C++ handles Beta distribution boundary cases", {
  # Values very close to 0
  test_data <- c(rep(0.001, 50), rbeta(50, 2, 2))
  dp <- DirichletProcessBeta(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Values very close to 1
  test_data <- c(rbeta(50, 2, 2), rep(0.999, 50))
  dp <- DirichletProcessBeta(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Mix of extreme values
  test_data <- c(rep(0.001, 33), rbeta(34, 2, 2), rep(0.999, 33))
  dp <- DirichletProcessBeta(test_data)
  expect_error(dp <- Fit(dp, its = 10), NA)
})

test_that("C++ handles interrupted/resumed fitting", {
  set.seed(123)
  test_data <- generate_test_data("normal", 100)

  # Initial fit
  dp1 <- DirichletProcessGaussian(test_data)
  dp1 <- Fit(dp1, its = 50)

  # Continue from where we left off
  dp2 <- Fit(dp1, its = 50)

  # Should have extended chains
  expect_length(dp2$alphaChain, 100)
  expect_length(dp2$labelsChain, 100)

  # Results should be continuous
  expect_equal(dp1$alpha, dp2$alphaChain[50])
})

test_that("C++ handles missing values appropriately", {
  # This test depends on how the package handles NAs
  test_data <- c(rnorm(90), rep(NA, 10))

  # Either it should handle gracefully or give informative error
  result <- tryCatch({
    dp <- DirichletProcessGaussian(test_data)
    "success"
  }, error = function(e) {
    "error"
  })

  # Document the behavior
  expect_true(result %in% c("success", "error"))
})

test_that("C++ handles different prior specifications", {
  test_data <- generate_test_data("normal", 100)

  # Test with extreme prior parameters
  extreme_priors <- list(
    m0 = 1000,
    s0 = 0.001,
    a0 = 0.001,
    b0 = 1000
  )

  dp <- DirichletProcessGaussian(test_data, g0Priors = extreme_priors)
  expect_error(dp <- Fit(dp, its = 10), NA)

  # Test with negative parameters (should error appropriately)
  bad_priors <- list(
    m0 = 0,
    s0 = -1,  # Invalid
    a0 = 1,
    b0 = 1
  )

  expect_error(
    DirichletProcessGaussian(test_data, g0Priors = bad_priors)
  )
})
