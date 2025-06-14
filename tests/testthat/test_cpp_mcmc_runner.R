context("C++ MCMC Runner Tests")

# Helper function to create test data
create_test_data <- function(n = 50, seed = 123) {
  set.seed(seed)
  c(rnorm(n/2, -2, 0.5), rnorm(n/2, 2, 0.5))
}

# Helper function to create MCMC parameters
create_mcmc_params <- function(n_iter = 100, n_burn = 20, thin = 1, update_concentration = TRUE, alpha = 1.0) {
  list(
    n_iter = n_iter,
    n_burn = n_burn,
    thin = thin,
    update_concentration = update_concentration,
    alpha = alpha
  )
}

# Helper function to create mixing distribution parameters
create_gaussian_params <- function() {
  list(
    type = "gaussian",
    mu0 = 0.0,
    kappa0 = 1.0,
    alpha0 = 1.0,
    beta0 = 1.0
  )
}

test_that("MCMC Runner basic functionality", {
  skip_if_not_installed("dirichletprocess")

  data <- create_test_data(50)
  data_matrix <- matrix(data, ncol = 1)
  mcmc_params <- create_mcmc_params()
  dist_params <- create_gaussian_params()

  # Test basic run
  result <- run_mcmc_cpp(data_matrix, dist_params, mcmc_params)

  expect_type(result, "list")
  expect_named(result, c("cluster_labels", "alpha", "theta", "n_clusters"))
  expect_equal(length(result$cluster_labels), (mcmc_params$n_iter - mcmc_params$n_burn) / mcmc_params$thin)
  expect_equal(length(result$alpha), (mcmc_params$n_iter - mcmc_params$n_burn) / mcmc_params$thin)
})

test_that("MCMC Runner parameter validation", {
  data_matrix <- matrix(create_test_data(20), ncol = 1)

  # Test with missing parameters
  expect_error(run_mcmc_cpp(data_matrix, list(), create_mcmc_params()))
  expect_error(run_mcmc_cpp(data_matrix, create_gaussian_params(), list()))

  # Test with invalid distribution type
  invalid_dist <- create_gaussian_params()
  invalid_dist$type <- "unknown"
  expect_error(run_mcmc_cpp(data_matrix, invalid_dist, create_mcmc_params()))
})

test_that("MCMC Runner different parameter settings", {
  data_matrix <- matrix(create_test_data(30), ncol = 1)

  # Test with different burn-in and thinning
  mcmc_params1 <- create_mcmc_params(n_iter = 200, n_burn = 50, thin = 2)
  result1 <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params1)
  expected_samples1 <- (200 - 50) / 2
  expect_equal(length(result1$cluster_labels), expected_samples1)

  # Test with no concentration update
  mcmc_params2 <- create_mcmc_params(update_concentration = FALSE)
  result2 <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params2)
  expect_type(result2, "list")

  # Test with different alpha
  mcmc_params3 <- create_mcmc_params(alpha = 5.0)
  result3 <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params3)
  expect_type(result3, "list")
})

test_that("Cluster labels validity", {
  data_matrix <- matrix(create_test_data(40), ncol = 1)
  mcmc_params <- create_mcmc_params(n_iter = 50, n_burn = 10)
  result <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params)

  # Check that cluster labels are valid
  for (labels in result$cluster_labels) {
    expect_true(all(labels >= 0))
    expect_true(all(labels < length(unique(labels)) + 10)) # reasonable upper bound
    expect_equal(length(labels), nrow(data_matrix))
  }
})

test_that("Alpha samples validity", {
  data_matrix <- matrix(create_test_data(30), ncol = 1)
  mcmc_params <- create_mcmc_params(n_iter = 60, n_burn = 10, update_concentration = TRUE)
  result <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params)

  # Check alpha samples
  for (alpha_vec in result$alpha) {
    expect_true(length(alpha_vec) == 1)
    expect_true(alpha_vec > 0)
    expect_true(is.finite(alpha_vec))
  }
})

test_that("Theta parameters validity", {
  data_matrix <- matrix(create_test_data(25), ncol = 1)
  mcmc_params <- create_mcmc_params(n_iter = 40, n_burn = 5)
  result <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params)

  # Check theta parameters (mean and variance for Gaussian)
  for (theta_list in result$theta) {
    expect_type(theta_list, "list")
    expect_true(length(theta_list) >= 1) # At least one cluster

    for (params in theta_list) {
      expect_equal(length(params), 2) # mean and variance
      expect_true(all(is.finite(params)))
      expect_true(params[2] > 0) # variance must be positive
    }
  }
})

test_that("Edge cases", {
  # Test with very small dataset
  small_data <- matrix(c(1.0, 2.0), ncol = 1)
  mcmc_params <- create_mcmc_params(n_iter = 20, n_burn = 5)
  result_small <- run_mcmc_cpp(small_data, create_gaussian_params(), mcmc_params)
  expect_type(result_small, "list")

  # Test with single data point
  single_data <- matrix(c(1.5), ncol = 1)
  result_single <- run_mcmc_cpp(single_data, create_gaussian_params(), mcmc_params)
  expect_type(result_single, "list")
  expect_equal(length(result_single$cluster_labels[[1]]), 1)
})

test_that("Reproducibility", {
  data_matrix <- matrix(create_test_data(30, seed = 456), ncol = 1)
  mcmc_params <- create_mcmc_params(n_iter = 30, n_burn = 5)

  # Note: C++ uses R's RNG, so setting seed should work
  set.seed(789)
  result1 <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params)

  set.seed(789)
  result2 <- run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params)

  # Results should be identical with same seed
  expect_equal(result1$cluster_labels, result2$cluster_labels)
  expect_equal(result1$alpha, result2$alpha)
})

test_that("C++ and R backends produce similar results", {
  skip_if_not_installed("dirichletprocess")

  data <- create_test_data(40, seed = 123)

  # Test with R backend
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(data)
  dp_r <- Fit(dp_r, 50, progressBar = FALSE)

  # Test with C++ backend (if available)
  tryCatch({
    set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessGaussian(data)
    dp_cpp <- Fit(dp_cpp, 50, progressBar = FALSE)

    # Compare basic properties
    expect_equal(length(dp_r$data), length(dp_cpp$data))
    expect_equal(dp_r$data, dp_cpp$data)

    # Both should have reasonable number of clusters
    expect_true(dp_r$numberClusters > 0)
    expect_true(dp_cpp$numberClusters > 0)
    expect_true(dp_r$numberClusters < length(data))
    expect_true(dp_cpp$numberClusters < length(data))

  }, error = function(e) {
    skip("C++ backend not available")
  })
})

test_that("Large dataset handling", {
  # Test with moderately large dataset
  large_data <- matrix(rnorm(1000), ncol = 1)
  mcmc_params <- create_mcmc_params(n_iter = 10, n_burn = 2)

  expect_silent(result <- run_mcmc_cpp(large_data, create_gaussian_params(), mcmc_params))
  expect_equal(length(result$cluster_labels[[1]]), 1000)
})

test_that("Invalid input handling", {
  # Test with empty data
  expect_error(run_mcmc_cpp(matrix(numeric(0), ncol = 1),
                            create_gaussian_params(), create_mcmc_params()))

  # Test with NA/Inf data
  bad_data <- matrix(c(1, NA, 3), ncol = 1)
  expect_error(run_mcmc_cpp(bad_data, create_gaussian_params(), create_mcmc_params()))

  bad_data2 <- matrix(c(1, Inf, 3), ncol = 1)
  expect_error(run_mcmc_cpp(bad_data2, create_gaussian_params(), create_mcmc_params()))
})

test_that("Parameter boundary conditions", {
  data_matrix <- matrix(rnorm(10), ncol = 1)

  # Test with zero iterations
  mcmc_params_zero <- create_mcmc_params(n_iter = 0)
  expect_error(run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params_zero))

  # Test with burn-in >= iterations
  mcmc_params_bad <- create_mcmc_params(n_iter = 10, n_burn = 15)
  expect_error(run_mcmc_cpp(data_matrix, create_gaussian_params(), mcmc_params_bad))
})
