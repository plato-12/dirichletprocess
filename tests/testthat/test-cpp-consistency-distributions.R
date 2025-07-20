# tests/testthat/test-cpp-consistency-distributions.R
#
# R/C++ Consistency Tests for All Distributions
# 
# This file tests statistical equivalence between R and C++ implementations
# across all 6 supported distributions with comprehensive validation.
#
# PERFORMANCE OPTIMIZATION:
# - DEV_MODE (default): Fast testing for development (50 iterations, smaller samples)
# - PRODUCTION_MODE: Full validation testing (200 iterations, larger samples)
#
# Usage:
# - Development: Sys.setenv(DP_DEV_TESTING = "TRUE") (default)
# - Production:  Sys.setenv(DP_DEV_TESTING = "FALSE")
#
# Speed improvement: ~75% faster in DEV_MODE while maintaining validation quality

# Development vs Production testing configuration
DEV_MODE <- Sys.getenv("DP_DEV_TESTING", "TRUE") == "TRUE"
BASE_ITERATIONS <- if (DEV_MODE) 50 else 200
BASE_SAMPLE_SIZE <- if (DEV_MODE) 50 else 100

test_that("Normal distribution R/C++ consistency", {
  # Generate test data
  set.seed(123)
  test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 0, 2), sd = 1)

  # Run consistency tests
  results <- validate_r_cpp_consistency("normal", test_data, iterations = BASE_ITERATIONS)

  # Assertions
  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
  expect_lt(results$param_max_diff, PARAM_TOLERANCE)
})

test_that("Exponential distribution R/C++ consistency", {
  set.seed(123)
  test_data <- rexp(BASE_SAMPLE_SIZE, rate = c(0.5, 1, 2))

  results <- validate_r_cpp_consistency("exponential", test_data, iterations = BASE_ITERATIONS)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("Beta distribution R/C++ consistency", {
  set.seed(123)
  beta_size <- if (DEV_MODE) 25 else 50
  test_data <- c(rbeta(beta_size, 2, 5), rbeta(beta_size, 5, 2))

  results <- validate_r_cpp_consistency("beta", test_data, iterations = BASE_ITERATIONS)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("Weibull distribution R/C++ consistency", {
  set.seed(123)
  test_data <- rweibull(BASE_SAMPLE_SIZE, shape = c(0.5, 1.5, 3), scale = 1)

  results <- validate_r_cpp_consistency("weibull", test_data, iterations = BASE_ITERATIONS)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("MVNormal distribution R/C++ consistency", {
  set.seed(123)
  mu1 <- c(0, 0)
  mu2 <- c(3, 3)
  sigma <- diag(2)
  mvn_size <- if (DEV_MODE) 25 else 50
  test_data <- rbind(
    mvtnorm::rmvnorm(mvn_size, mu1, sigma),
    mvtnorm::rmvnorm(mvn_size, mu2, sigma)
  )

  results <- validate_r_cpp_consistency("mvnormal", test_data, iterations = BASE_ITERATIONS)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("MVNormal2 distribution R/C++ consistency", {
  set.seed(123)
  mu1 <- c(-2, -2)
  mu2 <- c(2, 2)
  sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  mvn2_size <- if (DEV_MODE) 25 else 50
  test_data <- rbind(
    mvtnorm::rmvnorm(mvn2_size, mu1, sigma),
    mvtnorm::rmvnorm(mvn2_size, mu2, sigma)
  )

  results <- validate_r_cpp_consistency("mvnormal2", test_data, iterations = BASE_ITERATIONS)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

# Test different sample sizes
test_that("Consistency holds for different sample sizes", {
  sample_sizes <- if (DEV_MODE) c(25, 50) else c(50, 100, 500)
  dev_iterations <- if (DEV_MODE) 30 else 100
  dev_runs <- if (DEV_MODE) 2 else 3

  for (n in sample_sizes) {
    test_data <- generate_test_data("normal", n)
    results <- validate_r_cpp_consistency("normal", test_data, iterations = dev_iterations, n_runs = dev_runs)

    expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE,
              info = paste("Failed for sample size", n))
    expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE,
              info = paste("Failed for sample size", n))
  }
})

# Test different iteration counts
test_that("Consistency holds for different iteration counts", {
  test_data <- generate_test_data("exponential", BASE_SAMPLE_SIZE)
  iteration_counts <- if (DEV_MODE) c(25, 50) else c(50, 100, 500)
  iter_runs <- if (DEV_MODE) 2 else 3

  for (its in iteration_counts) {
    results <- validate_r_cpp_consistency("exponential", test_data,
                                          iterations = its, n_runs = iter_runs)

    expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE,
              info = paste("Failed for", its, "iterations"))
    expect_gt(results$likelihood_correlation, 0.9,  # Slightly lower for fewer iterations
              info = paste("Failed for", its, "iterations"))
  }
})
