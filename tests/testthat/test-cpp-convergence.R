# tests/testthat/test-cpp-convergence.R
#
# R/C++ Convergence Diagnostic Tests
# 
# This file tests convergence behavior and MCMC diagnostics between R and C++ implementations.
# These tests require longer chains and multiple runs for statistical validity.
#
# PERFORMANCE OPTIMIZATION:
# - DEV_MODE (default): Fast testing for development (reduced iterations, fewer chains)
# - PRODUCTION_MODE: Full convergence validation (long chains, multiple runs)
#
# Usage:
# - Development: Sys.setenv(DP_DEV_TESTING = "TRUE") (default)
# - Production:  Sys.setenv(DP_DEV_TESTING = "FALSE")
#
# Speed improvement: ~80% faster in DEV_MODE while maintaining convergence validation

library(coda)

# Development vs Production testing configuration
DEV_MODE <- Sys.getenv("DP_DEV_TESTING", "TRUE") == "TRUE"
CONV_ITERATIONS <- if (DEV_MODE) 250 else 1000
MULTI_ITERATIONS <- if (DEV_MODE) 150 else 500
CHAIN_COUNT <- if (DEV_MODE) 2 else 3
BASE_SAMPLE_SIZE <- if (DEV_MODE) 75 else 100
LARGE_SAMPLE_SIZE <- if (DEV_MODE) 100 else 200

test_that("R and C++ show similar convergence behavior", {
  set.seed(123)
  test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 2))

  # Longer chains for convergence analysis
  iterations <- CONV_ITERATIONS

  # R implementation
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(test_data)
  dp_r <- Fit(dp_r, its = iterations)

  # C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(test_data)
  dp_cpp <- Fit(dp_cpp, its = iterations)

  # Convert to mcmc objects for diagnostics
  mcmc_r <- mcmc(dp_r$alphaChain)
  mcmc_cpp <- mcmc(dp_cpp$alphaChain)

  # Effective sample size should be similar
  ess_r <- effectiveSize(mcmc_r)
  ess_cpp <- effectiveSize(mcmc_cpp)

  expect_equal(ess_r, ess_cpp, tolerance = 0.2)

  # Geweke diagnostics should both indicate convergence
  geweke_r <- geweke.diag(mcmc_r)$z
  geweke_cpp <- geweke.diag(mcmc_cpp)$z

  expect_lt(abs(geweke_r), 2)  # Within 2 standard deviations
  expect_lt(abs(geweke_cpp), 2)
})

test_that("Multiple chains show similar behavior", {
  set.seed(123)
  test_data <- generate_test_data("exponential", BASE_SAMPLE_SIZE)
  n_chains <- CHAIN_COUNT
  iterations <- MULTI_ITERATIONS

  # Run multiple chains for R
  r_chains <- list()
  for (i in 1:n_chains) {
    set.seed(123 + i)
    set_use_cpp(FALSE)
    dp <- DirichletProcessExponential(test_data)
    dp <- Fit(dp, its = iterations)
    r_chains[[i]] <- mcmc(dp$alphaChain)
  }
  r_mcmc_list <- mcmc.list(r_chains)

  # Run multiple chains for C++
  cpp_chains <- list()
  for (i in 1:n_chains) {
    set.seed(123 + i)
    set_use_cpp(TRUE)
    dp <- DirichletProcessExponential(test_data)
    dp <- Fit(dp, its = iterations)
    cpp_chains[[i]] <- mcmc(dp$alphaChain)
  }
  cpp_mcmc_list <- mcmc.list(cpp_chains)

  # Compare Gelman-Rubin statistics
  r_gelman <- gelman.diag(r_mcmc_list)
  cpp_gelman <- gelman.diag(cpp_mcmc_list)

  # Both should indicate convergence (close to 1)
  expect_lt(r_gelman$psrf[1], 1.1)
  expect_lt(cpp_gelman$psrf[1], 1.1)

  # Should be similar between R and C++
  expect_equal(r_gelman$psrf[1], cpp_gelman$psrf[1], tolerance = 0.1)
})

test_that("Autocorrelation patterns are similar", {
  set.seed(123)
  beta_size <- if (DEV_MODE) 75 else 150
  test_data <- generate_test_data("beta", beta_size)
  iterations <- CONV_ITERATIONS

  # R implementation
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessBeta(test_data)
  dp_r <- Fit(dp_r, its = iterations)

  # C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessBeta(test_data)
  dp_cpp <- Fit(dp_cpp, its = iterations)

  # Compute autocorrelations
  acf_r <- acf(dp_r$alphaChain, plot = FALSE)
  acf_cpp <- acf(dp_cpp$alphaChain, plot = FALSE)

  # Compare first 10 lags
  expect_equal(acf_r$acf[1:10], acf_cpp$acf[1:10], tolerance = 0.1)
})

test_that("Burn-in behavior is consistent", {
  set.seed(123)
  test_data <- generate_test_data("weibull", BASE_SAMPLE_SIZE)
  iterations <- MULTI_ITERATIONS

  # Track convergence metrics over time
  check_interval <- if (DEV_MODE) 30 else 50
  check_points <- seq(check_interval, iterations, by = check_interval)
  r_means <- numeric(length(check_points))
  cpp_means <- numeric(length(check_points))

  # R implementation
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessWeibull(test_data)

  for (i in seq_along(check_points)) {
    dp_r <- Fit(dp_r, its = ifelse(i == 1, check_points[1],
                                   check_points[i] - check_points[i-1]))
    r_means[i] <- mean(dp_r$alphaChain)
  }

  # C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessWeibull(test_data)

  for (i in seq_along(check_points)) {
    dp_cpp <- Fit(dp_cpp, its = ifelse(i == 1, check_points[1],
                                       check_points[i] - check_points[i-1]))
    cpp_means[i] <- mean(dp_cpp$alphaChain)
  }

  # Should converge to similar values
  expect_equal(r_means, cpp_means, tolerance = 0.1)

  # Both should stabilize (decreasing variance)
  r_diffs <- abs(diff(r_means))
  cpp_diffs <- abs(diff(cpp_means))

  expect_true(all(r_diffs[-1] <= r_diffs[-length(r_diffs)] + 0.01))
  expect_true(all(cpp_diffs[-1] <= cpp_diffs[-length(cpp_diffs)] + 0.01))
})

test_that("Posterior predictive distributions are similar", {
  set.seed(123)
  test_data <- generate_test_data("mvnormal", BASE_SAMPLE_SIZE)
  iterations <- MULTI_ITERATIONS
  n_posterior_samples <- if (DEV_MODE) 50 else 100

  # R implementation
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessMvnormal(test_data)
  dp_r <- Fit(dp_r, its = iterations)

  # C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessMvnormal(test_data)
  dp_cpp <- Fit(dp_cpp, its = iterations)

  # Generate posterior predictive samples
  set.seed(456)
  r_predictive <- PosteriorDraw(dp_r, n_posterior_samples)

  set.seed(456)
  cpp_predictive <- PosteriorDraw(dp_cpp, n_posterior_samples)

  # Compare distributions (using first dimension for simplicity)
  r_vals <- r_predictive[, 1]
  cpp_vals <- cpp_predictive[, 1]

  # Kolmogorov-Smirnov test - should not reject null hypothesis
  ks_test <- ks.test(r_vals, cpp_vals)
  expect_gt(ks_test$p.value, 0.05)

  # Compare summary statistics
  expect_equal(mean(r_vals), mean(cpp_vals), tolerance = 0.1)
  expect_equal(sd(r_vals), sd(cpp_vals), tolerance = 0.1)
})

test_that("Convergence diagnostics for cluster counts", {
  set.seed(123)
  test_data <- generate_test_data("normal", LARGE_SAMPLE_SIZE)
  iterations <- CONV_ITERATIONS

  # R implementation
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(test_data)
  dp_r <- Fit(dp_r, its = iterations)

  # C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(test_data)
  dp_cpp <- Fit(dp_cpp, its = iterations)

  # Extract cluster counts
  r_clusters <- sapply(dp_r$labelsChain, function(x) length(unique(x)))
  cpp_clusters <- sapply(dp_cpp$labelsChain, function(x) length(unique(x)))

  # Convert to mcmc objects
  mcmc_r_clusters <- mcmc(r_clusters)
  mcmc_cpp_clusters <- mcmc(cpp_clusters)

  # Both should show convergence
  geweke_r <- geweke.diag(mcmc_r_clusters)$z
  geweke_cpp <- geweke.diag(mcmc_cpp_clusters)$z

  expect_lt(abs(geweke_r), 2)
  expect_lt(abs(geweke_cpp), 2)

  # Should have similar posterior distributions
  expect_equal(mean(r_clusters), mean(cpp_clusters), tolerance = 0.5)
  expect_equal(sd(r_clusters), sd(cpp_clusters), tolerance = 0.2)
})
