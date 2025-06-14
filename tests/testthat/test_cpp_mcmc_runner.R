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

  # Only test if C++ implementation actually exists
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    skip("C++ MCMC implementation not compiled")
  }

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
  # Only test if C++ implementation actually exists
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    skip("C++ MCMC implementation not compiled")
  }

  data_matrix <- matrix(create_test_data(20), ncol = 1)

  # Test with missing parameters
  expect_error(run_mcmc_cpp(data_matrix, list(), create_mcmc_params()))
  expect_error(run_mcmc_cpp(data_matrix, create_gaussian_params(), list()))

  # Test with invalid distribution type
  invalid_dist <- create_gaussian_params()
  invalid_dist$type <- "unknown"
  expect_error(run_mcmc_cpp(data_matrix, invalid_dist, create_mcmc_params()))
})

test_that("R backend works correctly", {
  skip_if_not_installed("dirichletprocess")

  data <- create_test_data(40, seed = 123)

  # Test with R backend - this should ALWAYS work
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(data)

  # Use standard Fit method that works with R backend
  dp_r <- Fit(dp_r, 20, updatePrior = FALSE, progressBar = FALSE)

  # Verify R backend worked
  expect_true(exists("numberClusters", where = dp_r))
  expect_true(dp_r$numberClusters > 0)
  expect_equal(length(dp_r$data), length(data))
})

test_that("C++ backend switching", {
  skip_if_not_installed("dirichletprocess")

  # Only test C++ if it's actually available
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    skip("C++ MCMC implementation not compiled")
  }

  data <- create_test_data(30, seed = 123)

  # Test C++ backend
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(data)

  # This should either work with C++ or fall back to R
  expect_error({
    dp_cpp <- Fit(dp_cpp, 10, progressBar = FALSE)
  }, NA)  # Should not error

  expect_true(dp_cpp$numberClusters > 0)
})

test_that("Invalid input handling", {
  # Only test if C++ implementation actually exists
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    skip("C++ MCMC implementation not compiled")
  }

  # Test with empty data
  expect_error(run_mcmc_cpp(matrix(numeric(0), ncol = 1),
                            create_gaussian_params(), create_mcmc_params()))

  # Test with NA/Inf data
  bad_data <- matrix(c(1, NA, 3), ncol = 1)
  expect_error(run_mcmc_cpp(bad_data, create_gaussian_params(), create_mcmc_params()))

  bad_data2 <- matrix(c(1, Inf, 3), ncol = 1)
  expect_error(run_mcmc_cpp(bad_data2, create_gaussian_params(), create_mcmc_params()))
})

test_that("Backend switching works correctly", {
  skip_if_not_installed("dirichletprocess")

  data <- create_test_data(30, seed = 100)

  # Test R backend
  set_use_cpp(FALSE)
  expect_false(using_cpp())

  # Create DP object and fit with R backend
  dp_r <- DirichletProcessGaussian(data)
  expect_error({
    dp_r_fit <- Fit(dp_r, 10, updatePrior = FALSE, progressBar = FALSE)
  }, NA)  # Should not error

  # Test C++ backend status check
  cpp_status <- get_cpp_status()
  expect_type(cpp_status, "list")
  expect_true(all(sapply(cpp_status, is.logical)))
})

test_that("get_cpp_status function works", {
  # This should not error and should return a list
  status <- get_cpp_status()
  expect_type(status, "list")
  expect_true(all(sapply(status, is.logical)))

  # Should have expected components
  expected_components <- c("mcmc_runner", "gaussian_likelihood")
  expect_true(any(expected_components %in% names(status)))
})
