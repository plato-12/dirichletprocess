# tests/testthat/test-cpp-consistency-mvnormal.R
#
# R/C++ Consistency Tests for Multivariate Normal Distribution
# 
# This file tests statistical equivalence between R and C++ implementations
# for the Multivariate Normal distribution with comprehensive validation.

# Development vs Production testing configuration
DEV_MODE <- Sys.getenv("DP_DEV_TESTING", "TRUE") == "TRUE"
BASE_ITERATIONS <- if (DEV_MODE) 50 else 200

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