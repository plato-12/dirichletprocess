# tests/testthat/test-cpp-consistency-distributions.R

test_that("Normal distribution R/C++ consistency", {
  # Generate test data
  set.seed(123)
  test_data <- rnorm(100, mean = c(-2, 0, 2), sd = 1)

  # Run consistency tests
  results <- validate_r_cpp_consistency("normal", test_data, iterations = 200)

  # Assertions
  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
  expect_lt(results$param_max_diff, PARAM_TOLERANCE)
})

test_that("Exponential distribution R/C++ consistency", {
  set.seed(123)
  test_data <- rexp(100, rate = c(0.5, 1, 2))

  results <- validate_r_cpp_consistency("exponential", test_data, iterations = 200)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("Beta distribution R/C++ consistency", {
  set.seed(123)
  test_data <- c(rbeta(50, 2, 5), rbeta(50, 5, 2))

  results <- validate_r_cpp_consistency("beta", test_data, iterations = 200)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("Weibull distribution R/C++ consistency", {
  set.seed(123)
  test_data <- rweibull(100, shape = c(0.5, 1.5, 3), scale = 1)

  results <- validate_r_cpp_consistency("weibull", test_data, iterations = 200)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("MVNormal distribution R/C++ consistency", {
  set.seed(123)
  mu1 <- c(0, 0)
  mu2 <- c(3, 3)
  sigma <- diag(2)
  test_data <- rbind(
    mvtnorm::rmvnorm(50, mu1, sigma),
    mvtnorm::rmvnorm(50, mu2, sigma)
  )

  results <- validate_r_cpp_consistency("mvnormal", test_data, iterations = 200)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("MVNormal2 distribution R/C++ consistency", {
  set.seed(123)
  mu1 <- c(-2, -2)
  mu2 <- c(2, 2)
  sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  test_data <- rbind(
    mvtnorm::rmvnorm(50, mu1, sigma),
    mvtnorm::rmvnorm(50, mu2, sigma)
  )

  results <- validate_r_cpp_consistency("mvnormal2", test_data, iterations = 200)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

# Test different sample sizes
test_that("Consistency holds for different sample sizes", {
  sample_sizes <- c(50, 100, 500)

  for (n in sample_sizes) {
    test_data <- generate_test_data("normal", n)
    results <- validate_r_cpp_consistency("normal", test_data, iterations = 100, n_runs = 3)

    expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE,
              info = paste("Failed for sample size", n))
    expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE,
              info = paste("Failed for sample size", n))
  }
})

# Test different iteration counts
test_that("Consistency holds for different iteration counts", {
  test_data <- generate_test_data("exponential", 100)
  iteration_counts <- c(50, 100, 500)

  for (its in iteration_counts) {
    results <- validate_r_cpp_consistency("exponential", test_data,
                                          iterations = its, n_runs = 3)

    expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE,
              info = paste("Failed for", its, "iterations"))
    expect_gt(results$likelihood_correlation, 0.9,  # Slightly lower for fewer iterations
              info = paste("Failed for", its, "iterations"))
  }
})
