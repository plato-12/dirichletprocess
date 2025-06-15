# Comprehensive unit tests for Gaussian MCMC C++ implementation
library(testthat)
library(dirichletprocess)

context("Gaussian MCMC C++ Implementation")

# Helper function to skip tests if C++ not available
skip_if_no_cpp <- function() {
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    skip("C++ MCMC implementation not compiled")
  }
}

# Helper function to generate test data
generate_gaussian_mixture <- function(n = 100, k = 3, dim = 1, seed = 123) {
  set.seed(seed)

  if (k == 2 && dim == 1) {
    # Simple two-component mixture for most tests
    data <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 2, 0.5))
  } else if (k == 3 && dim == 1) {
    # Three-component mixture
    data <- c(rnorm(n/3, -3, 0.5), rnorm(n/3, 0, 0.5), rnorm(n/3, 3, 0.5))
  } else {
    # General case
    components <- sample(1:k, n, replace = TRUE)
    means <- seq(-3, 3, length.out = k)
    data <- numeric(n)
    for (i in 1:n) {
      data[i] <- rnorm(1, means[components[i]], 0.5)
    }
  }

  return(data)
}

test_that("C++ backend can be enabled and disabled", {
  # Test enabling C++ backend
  set_use_cpp(TRUE)
  expect_true(using_cpp())

  # Test disabling C++ backend
  set_use_cpp(FALSE)
  expect_false(using_cpp())

  # Check status function
  status <- get_cpp_status()
  expect_type(status, "list")
  expect_true(all(sapply(status, is.logical)))
})

test_that("DirichletProcessGaussian initializes correctly", {
  test_data <- generate_gaussian_mixture(n = 50, k = 2)

  # Test with default priors
  dp <- DirichletProcessGaussian(test_data)

  expect_s3_class(dp, "dirichletprocess")
  expect_equal(length(dp$data), 50)
  expect_equal(dp$n, 50)
  expect_true(!is.null(dp$mixingDistribution))
  expect_true(inherits(dp$mixingDistribution, "normal"))

  # Test with custom priors
  dp_custom <- DirichletProcessGaussian(
    test_data,
    g0Priors = c(0, 2, 2, 2),
    alphaPriors = c(1, 1)
  )

  expect_equal(dp_custom$mixingDistribution$priorParameters, c(0, 2, 2, 2))
  expect_equal(dp_custom$alphaPriorParameters, c(1, 1))
})

test_that("MCMC runner works with C++ backend", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 100, k = 2)

  # Enable C++ backend
  set_use_cpp(TRUE)

  # Create and fit model - using 'its' parameter name
  dp <- DirichletProcessGaussian(test_data)
  dp_fit <- Fit(dp, its = 100, progressBar = FALSE)

  # Check basic structure
  expect_true(!is.null(dp_fit$clusterLabels))
  expect_true(!is.null(dp_fit$clusterParameters))
  expect_true(!is.null(dp_fit$numberClusters))
  expect_equal(length(dp_fit$clusterLabels), 100)

  # Check that clusters were found
  expect_true(dp_fit$numberClusters >= 1)
  expect_true(dp_fit$numberClusters <= 100)  # Can't have more clusters than data points
})

test_that("C++ and R backends produce comparable results", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 50, k = 2, seed = 456)

  # Run with R backend
  set_use_cpp(FALSE)
  set.seed(789)
  dp_r <- DirichletProcessGaussian(test_data)
  dp_r_fit <- Fit(dp_r, its = 50, progressBar = FALSE)

  # Run with C++ backend
  set_use_cpp(TRUE)
  set.seed(789)
  dp_cpp <- DirichletProcessGaussian(test_data)
  dp_cpp_fit <- Fit(dp_cpp, its = 50, progressBar = FALSE)

  # Results won't be identical due to implementation differences,
  # but should find similar number of clusters
  expect_true(abs(dp_r_fit$numberClusters - dp_cpp_fit$numberClusters) <= 2)

  # Both should identify the bimodal structure (2-4 clusters typically)
  expect_true(dp_r_fit$numberClusters >= 1 && dp_r_fit$numberClusters <= 5)
  expect_true(dp_cpp_fit$numberClusters >= 1 && dp_cpp_fit$numberClusters <= 5)
})

test_that("Cluster parameters are updated correctly", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 80, k = 2)

  set_use_cpp(TRUE)
  dp <- DirichletProcessGaussian(test_data)
  dp_fit <- Fit(dp, its = 100, progressBar = FALSE)

  # Check cluster parameters structure
  params <- dp_fit$clusterParameters
  n_clusters <- dp_fit$numberClusters

  # Should have parameters for each cluster
  expect_type(params, "list")
  expect_equal(length(params), n_clusters)

  # Each cluster should have mean and variance
  for (i in 1:n_clusters) {
    expect_true(!is.null(params[[i]]$mu))
    expect_true(!is.null(params[[i]]$sig))
    expect_true(is.numeric(params[[i]]$mu))
    expect_true(is.numeric(params[[i]]$sig))
    expect_true(params[[i]]$sig > 0)  # Variance must be positive
  }
})

test_that("Alpha (concentration parameter) updates work", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 100, k = 3)

  set_use_cpp(TRUE)

  # Test with alpha updates enabled (default)
  dp1 <- DirichletProcessGaussian(test_data, alphaPriors = c(2, 4))
  dp1_fit <- Fit(dp1, its = 100, progressBar = FALSE)

  # Alpha should be positive and reasonable
  expect_true(dp1_fit$alpha > 0)
  expect_true(dp1_fit$alpha < 100)  # Sanity check

  # Check alpha chain if available
  if (!is.null(dp1_fit$alphaChain)) {
    expect_equal(length(dp1_fit$alphaChain), 100)
    expect_true(all(dp1_fit$alphaChain > 0))
  }
})

test_that("MCMC handles edge cases correctly", {
  skip_if_no_cpp()

  set_use_cpp(TRUE)

  # Test 1: Very small dataset
  small_data <- rnorm(5)
  dp_small <- DirichletProcessGaussian(small_data)
  expect_error(dp_small_fit <- Fit(dp_small, its = 20, progressBar = FALSE), NA)
  expect_true(dp_small_fit$numberClusters >= 1)
  expect_true(dp_small_fit$numberClusters <= 5)

  # Test 2: Identical data points
  identical_data <- rep(1.0, 20)
  dp_identical <- DirichletProcessGaussian(identical_data)
  expect_error(dp_identical_fit <- Fit(dp_identical, its = 20, progressBar = FALSE), NA)
  # Should typically find 1 cluster for identical data
  expect_true(dp_identical_fit$numberClusters >= 1)
  expect_true(dp_identical_fit$numberClusters <= 3)

  # Test 3: Data with outliers
  outlier_data <- c(rnorm(45, 0, 1), c(-10, 10, -10, 10, 15))
  dp_outlier <- DirichletProcessGaussian(outlier_data)
  expect_error(dp_outlier_fit <- Fit(dp_outlier, its = 50, progressBar = FALSE), NA)
})

test_that("Parameter validation works correctly", {
  skip_if_no_cpp()

  # Only test direct C++ function if it exists
  if (exists("run_mcmc_cpp")) {
    # Test invalid data inputs
    expect_error(run_mcmc_cpp(matrix(numeric(0), ncol = 1),
                              list(type = "gaussian", mu0 = 0, kappa0 = 1, alpha0 = 1, beta0 = 1),
                              list(n_iter = 10, n_burn = 0, thin = 1, update_concentration = TRUE, alpha = 1)))

    # Test with NA data
    bad_data <- matrix(c(1, NA, 3), ncol = 1)
    expect_error(run_mcmc_cpp(bad_data,
                              list(type = "gaussian", mu0 = 0, kappa0 = 1, alpha0 = 1, beta0 = 1),
                              list(n_iter = 10, n_burn = 0, thin = 1, update_concentration = TRUE, alpha = 1)))
  }
})

test_that("MCMC parameters are respected", {
  skip_if_no_cpp()

  # Only test if direct C++ function exists
  if (!exists("run_mcmc_cpp")) {
    skip("Direct C++ interface not available")
  }

  test_data <- generate_gaussian_mixture(n = 50, k = 2)
  data_matrix <- matrix(test_data, ncol = 1)

  # Test with different MCMC settings
  mixing_params <- list(type = "gaussian", mu0 = 0, kappa0 = 1, alpha0 = 1, beta0 = 1)

  # Test burn-in
  result1 <- run_mcmc_cpp(data_matrix, mixing_params,
                          list(n_iter = 100, n_burn = 20, thin = 1,
                               update_concentration = TRUE, alpha = 1))

  # Check if results have expected structure
  if (!is.null(result1$cluster_labels) && length(result1$cluster_labels) > 0) {
    expect_equal(length(result1$cluster_labels), 80)  # 100 - 20 burn-in
  }

  # Test thinning
  result2 <- run_mcmc_cpp(data_matrix, mixing_params,
                          list(n_iter = 100, n_burn = 0, thin = 5,
                               update_concentration = TRUE, alpha = 1))

  if (!is.null(result2$cluster_labels) && length(result2$cluster_labels) > 0) {
    expect_equal(length(result2$cluster_labels), 20)  # 100 / 5 thinning
  }
})

test_that("Likelihood calculations are correct", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 100, k = 2)

  set_use_cpp(TRUE)
  dp <- DirichletProcessGaussian(test_data)
  dp_fit <- Fit(dp, its = 50, progressBar = FALSE)

  # Check if likelihood values are stored
  if (!is.null(dp_fit$likelihoodChain)) {
    # Likelihood should be finite and negative (log-likelihood)
    expect_true(all(is.finite(dp_fit$likelihoodChain)))
    expect_true(all(dp_fit$likelihoodChain < 0))

    # Likelihood should generally improve over iterations
    first_quarter <- mean(dp_fit$likelihoodChain[1:12])
    last_quarter <- mean(dp_fit$likelihoodChain[38:50])
    # Allow for some stochasticity
    expect_true(last_quarter >= first_quarter * 0.95)
  }
})

test_that("Predictive distribution works with C++ backend", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 80, k = 2)

  set_use_cpp(TRUE)
  dp <- DirichletProcessGaussian(test_data)
  dp_fit <- Fit(dp, its = 100, progressBar = FALSE)

  # Test posterior clusters on new data
  new_data <- matrix(c(-2, 0, 2), ncol = 1)

  # This should work regardless of backend
  expect_error({
    pred_clusters <- PosteriorClusters(dp_fit, new_data)
  }, NA)

  if (exists("pred_clusters")) {
    expect_equal(nrow(pred_clusters), 3)
    expect_true(all(pred_clusters >= 0))
    expect_true(all(is.finite(pred_clusters)))
  }
})

test_that("Memory usage is reasonable", {
  skip_if_no_cpp()
  skip_on_cran()  # Memory tests can be flaky on CRAN

  test_data <- generate_gaussian_mixture(n = 500, k = 3)

  set_use_cpp(TRUE)

  # Get initial memory
  gc()
  mem_before <- gc()[2, 2]  # Used memory in MB

  # Run MCMC
  dp <- DirichletProcessGaussian(test_data)
  dp_fit <- Fit(dp, its = 100, progressBar = FALSE)

  # Check memory after
  gc()
  mem_after <- gc()[2, 2]

  # Memory increase should be reasonable (less than 50 MB for this test)
  mem_increase <- mem_after - mem_before
  expect_true(mem_increase < 50)
})

test_that("Prior specifications are respected", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 60, k = 2)

  set_use_cpp(TRUE)

  # Test with informative prior centered at 0
  dp_informative <- DirichletProcessGaussian(
    test_data,
    g0Priors = c(0, 10, 5, 5),  # Strong prior on mean near 0
    alphaPriors = c(1, 1)
  )
  dp_informative_fit <- Fit(dp_informative, its = 50, progressBar = FALSE)

  # Test with weak prior
  dp_weak <- DirichletProcessGaussian(
    test_data,
    g0Priors = c(0, 0.01, 1, 1),  # Weak prior on mean
    alphaPriors = c(1, 1)
  )
  dp_weak_fit <- Fit(dp_weak, its = 50, progressBar = FALSE)

  # Both should find clusters, but the number might differ
  expect_true(dp_informative_fit$numberClusters >= 1)
  expect_true(dp_weak_fit$numberClusters >= 1)
})

test_that("Cluster label updates work correctly", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 30, k = 2)

  set_use_cpp(TRUE)
  dp <- DirichletProcessGaussian(test_data)

  # Store initial state
  initial_labels <- dp$clusterLabels

  # Fit for a few iterations
  dp_fit <- Fit(dp, its = 10, progressBar = FALSE)

  # Labels should have changed
  final_labels <- dp_fit$clusterLabels
  expect_false(all(initial_labels == final_labels))

  # All labels should be valid
  expect_true(all(final_labels >= 1))
  expect_true(all(final_labels <= 30))
  expect_equal(length(final_labels), 30)
})

test_that("Consistency between multiple runs with same seed", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 40, k = 2)

  set_use_cpp(TRUE)

  # Run 1
  set.seed(12345)
  dp1 <- DirichletProcessGaussian(test_data)
  dp1_fit <- Fit(dp1, its = 30, progressBar = FALSE)

  # Run 2 with same seed
  set.seed(12345)
  dp2 <- DirichletProcessGaussian(test_data)
  dp2_fit <- Fit(dp2, its = 30, progressBar = FALSE)

  # Key results should be identical
  expect_equal(dp1_fit$numberClusters, dp2_fit$numberClusters)
  expect_equal(dp1_fit$alpha, dp2_fit$alpha)
  expect_equal(dp1_fit$clusterLabels, dp2_fit$clusterLabels)
})

test_that("Integration with other package functions works", {
  skip_if_no_cpp()

  test_data <- generate_gaussian_mixture(n = 100, k = 2)

  set_use_cpp(TRUE)
  dp <- DirichletProcessGaussian(test_data)
  dp_fit <- Fit(dp, its = 50, progressBar = FALSE)

  # Test that standard S3 methods work
  expect_error({
    print(dp_fit)
  }, NA)

  # Test plotting (if available)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_error({
      p <- plot(dp_fit)
    }, NA)
  }

  # Test cluster change functions
  expect_error({
    dp_changed <- ClusterLabelChange(dp_fit, 1, 2, 1)
  }, NA)
})
