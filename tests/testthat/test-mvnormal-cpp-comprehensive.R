# test_mvnormal_cpp_comprehensive.R
# Comprehensive unit tests for MVNormal C++ implementations

library(testthat)
library(dirichletprocess)
library(mvtnorm)

context("MVNormal Distribution C++ Implementation Comprehensive Tests")

# ==============================================================================
# Test Setup and Helper Functions
# ==============================================================================

create_test_data <- function(n = 100, d = 2, n_clusters = 2, seed = 123) {
  set.seed(seed)
  if (n_clusters == 1) {
    mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = diag(d))
  } else if (n_clusters == 2) {
    rbind(
      mvtnorm::rmvnorm(n/2, mean = rep(-2, d), sigma = diag(d) * 0.5),
      mvtnorm::rmvnorm(n/2, mean = rep(2, d), sigma = diag(d) * 1.5)
    )
  } else {
    stop("Only 1 or 2 clusters supported in test data generator")
  }
}

create_test_priors <- function(d = 2) {
  list(
    mu0 = rep(0, d),
    Lambda = diag(d),
    kappa0 = 1,
    nu = d + 2
  )
}

# ==============================================================================
# Test 1: Basic C++ Function Availability
# ==============================================================================

test_that("MVNormal C++ functions are available", {
  # Check core MCMC runner
  expect_true(exists("_dirichletprocess_run_mcmc_cpp"),
              "C++ MCMC runner not found")

  # Check MVNormal-specific conjugate functions
  expect_true(exists("conjugate_mvnormal_cluster_component_update_cpp"),
              "MVNormal cluster component update not found")
  expect_true(exists("conjugate_mvnormal_cluster_parameter_update_cpp"),
              "MVNormal cluster parameter update not found")

  # Check MVNormal distribution functions
  expect_true(exists("mvnormal_prior_draw_cpp"),
              "MVNormal prior draw function not found")
  expect_true(exists("mvnormal_posterior_draw_cpp"),
              "MVNormal posterior draw function not found")
  expect_true(exists("mvnormal_posterior_parameters_cpp"),
              "MVNormal posterior parameters function not found")
  expect_true(exists("mvnormal_predictive_cpp"),
              "MVNormal predictive function not found")
  expect_true(exists("mvnormal_likelihood_cpp"),
              "MVNormal likelihood function not found")
})

# ==============================================================================
# Test 2: Prior Distribution Tests
# ==============================================================================

test_that("MVNormal prior draw works correctly", {
  skip_if_not(exists("mvnormal_prior_draw_cpp"))

  # Test 2D case
  priors_2d <- create_test_priors(d = 2)

  # Single draw
  set.seed(456)
  draw_single <- mvnormal_prior_draw_cpp(priors_2d, n = 1)

  expect_type(draw_single, "list")
  expect_equal(length(draw_single$mu), 2)  # 2D mean
  expect_equal(dim(draw_single$sig), c(2, 2))  # 2x2 precision matrix

  # Multiple draws
  set.seed(789)
  draws_multiple <- mvnormal_prior_draw_cpp(priors_2d, n = 10)

  expect_equal(length(draws_multiple$mu), 20)  # 10 draws * 2 dimensions
  expect_equal(dim(draws_multiple$sig), c(2, 2, 10))  # 2x2x10 array

  # Test different dimensions
  for (d in c(1, 3, 5)) {
    priors_d <- create_test_priors(d = d)
    draw_d <- mvnormal_prior_draw_cpp(priors_d, n = 1)

    expect_equal(length(draw_d$mu), d)
    expect_equal(dim(draw_d$sig), c(d, d))
  }

  # Verify draws are from Wishart distribution (precision matrices should be positive definite)
  for (i in 1:10) {
    sig_i <- draws_multiple$sig[,,i]
    eigenvals <- eigen(sig_i, only.values = TRUE)$values
    expect_true(all(eigenvals > 0),
                info = sprintf("Draw %d: precision matrix not positive definite", i))
  }
})

# ==============================================================================
# Test 3: Posterior Distribution Tests
# ==============================================================================

test_that("MVNormal posterior parameters are computed correctly", {
  skip_if_not(exists("mvnormal_posterior_parameters_cpp"))

  # Test with simple data
  priors <- create_test_priors(d = 2)
  data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, byrow = TRUE)

  post_params <- mvnormal_posterior_parameters_cpp(priors, data)

  expect_type(post_params, "list")
  expect_true(all(c("mu_n", "kappa_n", "nu_n", "t_n") %in% names(post_params)))

  # Check dimensions
  expect_equal(length(post_params$mu_n), 2)
  expect_equal(dim(post_params$t_n), c(2, 2))
  expect_type(post_params$kappa_n, "double")
  expect_type(post_params$nu_n, "double")

  # Verify parameter updates follow conjugate update rules
  n_data <- nrow(data)
  expect_equal(post_params$kappa_n, priors$kappa0 + n_data)
  expect_equal(post_params$nu_n, priors$nu + n_data)

  # Test with edge cases
  # Single observation
  data_single <- matrix(c(0, 0), ncol = 2)
  post_single <- mvnormal_posterior_parameters_cpp(priors, data_single)
  expect_equal(post_single$kappa_n, priors$kappa0 + 1)

  # Zero-centered data
  data_zero <- matrix(rnorm(20), ncol = 2)
  data_zero <- scale(data_zero, center = TRUE, scale = FALSE)
  post_zero <- mvnormal_posterior_parameters_cpp(priors, data_zero)
  expect_true(all(abs(colMeans(data_zero)) < 1e-10))
})

test_that("MVNormal posterior draw works correctly", {
  skip_if_not(exists("mvnormal_posterior_draw_cpp"))

  priors <- create_test_priors(d = 2)
  data <- create_test_data(n = 30, d = 2, n_clusters = 1, seed = 111)

  # Single draw
  set.seed(222)
  post_draw <- mvnormal_posterior_draw_cpp(priors, data, n = 1)

  expect_type(post_draw, "list")
  expect_equal(length(post_draw$mu), 2)
  expect_equal(dim(post_draw$sig), c(2, 2))

  # Multiple draws
  set.seed(333)
  post_draws <- mvnormal_posterior_draw_cpp(priors, data, n = 100)

  expect_equal(length(post_draws$mu), 200)  # 100 draws * 2 dimensions
  expect_equal(dim(post_draws$sig), c(2, 2, 100))

  # Check that posterior means are reasonable (should be near data mean for large n)
  data_mean <- colMeans(data)
  posterior_means <- matrix(post_draws$mu, ncol = 2, byrow = TRUE)
  mean_of_posterior_means <- colMeans(posterior_means)

  # With 30 observations and reasonable prior, posterior mean should be close to data mean
  expect_true(all(abs(mean_of_posterior_means - data_mean) < 0.5),
              info = "Posterior mean too far from data mean")
})

# ==============================================================================
# Test 4: Likelihood Function Tests
# ==============================================================================

test_that("MVNormal likelihood function works correctly", {
  skip_if_not(exists("mvnormal_likelihood_cpp"))

  # Test 2D case
  mu <- c(0, 0)
  sigma <- diag(2)
  x <- matrix(c(0, 0, 1, 1, -1, -1), ncol = 2, byrow = TRUE)

  lik <- mvnormal_likelihood_cpp(x, mu, sigma)

  expect_type(lik, "double")
  expect_equal(length(lik), 3)
  expect_true(all(lik > 0))
  expect_true(all(lik <= 1/(2*pi)))  # Maximum density for 2D standard normal

  # Compare with mvtnorm package
  lik_mvtnorm <- dmvnorm(x, mean = mu, sigma = sigma)
  expect_equal(lik, lik_mvtnorm, tolerance = 1e-10)

  # Test with different covariance structures
  # Diagonal covariance
  sigma_diag <- diag(c(2, 0.5))
  lik_diag <- mvnormal_likelihood_cpp(x, mu, sigma_diag)
  lik_diag_ref <- dmvnorm(x, mean = mu, sigma = sigma_diag)
  expect_equal(lik_diag, lik_diag_ref, tolerance = 1e-10)

  # Full covariance
  sigma_full <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  lik_full <- mvnormal_likelihood_cpp(x, mu, sigma_full)
  lik_full_ref <- dmvnorm(x, mean = mu, sigma = sigma_full)
  expect_equal(lik_full, lik_full_ref, tolerance = 1e-10)

  # Test edge cases
  # Single point
  x_single <- matrix(c(0, 0), ncol = 2)
  lik_single <- mvnormal_likelihood_cpp(x_single, mu, sigma)
  expect_equal(length(lik_single), 1)

  # Extreme values (should have very low likelihood)
  x_extreme <- matrix(c(10, 10), ncol = 2)
  lik_extreme <- mvnormal_likelihood_cpp(x_extreme, mu, sigma)
  expect_true(lik_extreme < 1e-10)
})

# ==============================================================================
# Test 5: Predictive Distribution Tests
# ==============================================================================

test_that("MVNormal predictive distribution works correctly", {
  skip_if_not(exists("mvnormal_predictive_cpp"))

  priors <- create_test_priors(d = 2)
  x_new <- matrix(c(0, 0, 1, 1, -1, -1), ncol = 2, byrow = TRUE)

  pred <- mvnormal_predictive_cpp(priors, x_new)

  expect_type(pred, "double")
  expect_equal(length(pred), 3)
  expect_true(all(pred > 0))

  # Predictive should integrate likelihood over parameter uncertainty
  # For conjugate prior, this gives a multivariate t-distribution
  # Values at prior mean should have highest predictive probability
  x_at_mean <- matrix(priors$mu0, ncol = 2)
  pred_at_mean <- mvnormal_predictive_cpp(priors, x_at_mean)

  x_far <- matrix(c(5, 5), ncol = 2)
  pred_far <- mvnormal_predictive_cpp(priors, x_far)

  expect_true(pred_at_mean > pred_far)
})

# ==============================================================================
# Test 6: Cluster Update Functions (Conjugate)
# ==============================================================================

test_that("MVNormal cluster component update works correctly", {
  skip_if_not(exists("conjugate_mvnormal_cluster_component_update_cpp"))

  # Create a simple DP object
  data <- create_test_data(n = 50, d = 2, n_clusters = 2, seed = 444)
  priors <- create_test_priors(d = 2)
  dp <- DirichletProcessMvnormal(data, priors)
  dp <- Initialise(dp, numInitialClusters = 3)

  # Ensure arrays have sufficient size
  if (dim(dp$clusterParameters$mu)[3] < 10) {
    new_mu <- array(NA_real_, dim = c(1, 2, 10))
    new_sig <- array(NA_real_, dim = c(2, 2, 10))

    old_dim <- dim(dp$clusterParameters$mu)[3]
    new_mu[, , 1:old_dim] <- dp$clusterParameters$mu
    new_sig[, , 1:old_dim] <- dp$clusterParameters$sig

    dp$clusterParameters$mu <- new_mu
    dp$clusterParameters$sig <- new_sig
  }

  # Convert to 0-indexed for C++
  dp$clusterLabels <- dp$clusterLabels - 1

  # Run cluster component update
  result <- conjugate_mvnormal_cluster_component_update_cpp(dp)

  expect_type(result, "list")
  expect_true(all(c("clusterLabels", "pointsPerCluster", "numberClusters",
                    "clusterParameters") %in% names(result)))

  # Convert back to 1-indexed
  result$clusterLabels <- result$clusterLabels + 1

  # Validity checks
  expect_equal(length(result$clusterLabels), nrow(data))
  expect_true(all(result$clusterLabels >= 1))
  expect_true(all(result$clusterLabels <= result$numberClusters))
  expect_equal(sum(result$pointsPerCluster), nrow(data))
  expect_true(result$numberClusters >= 1)
  expect_true(result$numberClusters <= nrow(data))
})

test_that("MVNormal cluster parameter update works correctly", {
  skip_if_not(exists("conjugate_mvnormal_cluster_parameter_update_cpp"))

  # Create DP object with known cluster assignments
  data <- create_test_data(n = 40, d = 2, n_clusters = 2, seed = 555)
  priors <- create_test_priors(d = 2)
  dp <- DirichletProcessMvnormal(data, priors)
  dp <- Initialise(dp, numInitialClusters = 2)

  # Set cluster labels manually (first 20 in cluster 1, next 20 in cluster 2)
  dp$clusterLabels <- c(rep(1, 20), rep(2, 20))
  dp$pointsPerCluster <- c(20, 20)
  dp$numberClusters <- 2

  # Convert to 0-indexed
  dp$clusterLabels <- dp$clusterLabels - 1

  # Run parameter update
  updated_params <- conjugate_mvnormal_cluster_parameter_update_cpp(dp)

  expect_type(updated_params, "list")
  expect_true(all(c("mu", "sig") %in% names(updated_params)))

  # Check dimensions
  expect_equal(dim(updated_params$mu), c(1, 2, dim(dp$clusterParameters$mu)[3]))
  expect_equal(dim(updated_params$sig), c(2, 2, dim(dp$clusterParameters$sig)[3]))

  # Parameters for active clusters should be updated
  # Check that first two clusters have reasonable values
  mu1 <- updated_params$mu[1, , 1]
  mu2 <- updated_params$mu[1, , 2]

  # Means should be different for the two clusters
  expect_true(sqrt(sum((mu1 - mu2)^2)) > 0.1)
})

# ==============================================================================
# Test 7: Full MCMC Integration Tests
# ==============================================================================

test_that("MVNormal C++ MCMC produces statistically valid results", {
  skip_if_not(can_use_cpp(DirichletProcessMvnormal(matrix(1), list())))

  # Generate data with known structure
  set.seed(666)
  true_means <- list(c(-3, -3), c(0, 0), c(3, 3))
  true_covs <- list(diag(2) * 0.5, diag(2), diag(2) * 0.5)

  data <- rbind(
    mvtnorm::rmvnorm(30, true_means[[1]], true_covs[[1]]),
    mvtnorm::rmvnorm(40, true_means[[2]], true_covs[[2]]),
    mvtnorm::rmvnorm(30, true_means[[3]], true_covs[[3]])
  )

  priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2) * 2,
    kappa0 = 1,
    nu = 4
  )

  # Fit with C++ backend
  set_use_cpp(TRUE)
  set.seed(777)
  dp <- DirichletProcessMvnormal(data, priors)
  dp <- Fit(dp, 200, progressBar = FALSE)

  # Should find approximately 3 clusters (allowing for MCMC variability)
  expect_true(dp$numberClusters >= 2)
  expect_true(dp$numberClusters <= 5)

  # Check convergence indicators
  # Alpha chain should stabilize
  alpha_second_half <- dp$alphaChain[101:200]
  expect_true(sd(alpha_second_half) < 2)

  # Likelihood should generally increase and stabilize
  lik_second_half <- dp$likelihoodChain[101:200]
  expect_true(mean(lik_second_half) > mean(dp$likelihoodChain[1:50]))
})

# ==============================================================================
# Test 8: Edge Cases and Error Handling
# ==============================================================================

test_that("MVNormal C++ handles edge cases gracefully", {
  # 1D data (degenerate multivariate case)
  data_1d <- matrix(rnorm(20), ncol = 1)
  priors_1d <- list(
    mu0 = 0,
    Lambda = matrix(1, 1, 1),
    kappa0 = 1,
    nu = 2
  )

  # Should work without errors
  pred_1d <- mvnormal_predictive_cpp(priors_1d, data_1d)
  expect_equal(length(pred_1d), 20)

  # Small sample size
  data_small <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  priors_2d <- create_test_priors(d = 2)

  post_params_small <- mvnormal_posterior_parameters_cpp(priors_2d, data_small)
  expect_equal(post_params_small$nu_n, priors_2d$nu + 2)

  # Extreme correlations
  sigma_extreme <- matrix(c(1, 0.99, 0.99, 1), 2, 2)
  x_test <- matrix(c(0, 0), ncol = 2)

  # Should handle near-singular matrices
  lik_extreme <- mvnormal_likelihood_cpp(x_test, c(0, 0), sigma_extreme)
  expect_true(is.finite(lik_extreme))

  # Large dimension
  d_large <- 10
  priors_large <- create_test_priors(d = d_large)
  data_large <- matrix(rnorm(50 * d_large), ncol = d_large)

  # Should handle higher dimensions
  post_large <- mvnormal_posterior_draw_cpp(priors_large, data_large, n = 1)
  expect_equal(length(post_large$mu), d_large)
  expect_equal(dim(post_large$sig), c(d_large, d_large))
})

# ==============================================================================
# Test 9: Performance Benchmarks
# ==============================================================================

test_that("MVNormal C++ is faster than R implementation", {
  skip_if_not(can_use_cpp(DirichletProcessMvnormal(matrix(1), list())))
  skip_on_cran()  # Skip performance tests on CRAN

  # Create moderate-sized dataset
  data <- create_test_data(n = 200, d = 3, n_clusters = 2, seed = 888)
  priors <- create_test_priors(d = 3)

  # Time R implementation
  set_use_cpp(FALSE)
  time_r <- system.time({
    dp_r <- DirichletProcessMvnormal(data, priors)
    dp_r <- Fit(dp_r, 50, progressBar = FALSE)
  })["elapsed"]

  # Time C++ implementation
  set_use_cpp(TRUE)
  time_cpp <- system.time({
    dp_cpp <- DirichletProcessMvnormal(data, priors)
    dp_cpp <- Fit(dp_cpp, 50, progressBar = FALSE)
  })["elapsed"]

  speedup <- time_r / time_cpp

  cat("\nMVNormal Performance:")
  cat("\n  R time:", round(time_r, 3), "seconds")
  cat("\n  C++ time:", round(time_cpp, 3), "seconds")
  cat("\n  Speedup:", round(speedup, 1), "x\n")

  # Expect at least 5x speedup
  expect_true(speedup > 5,
              info = sprintf("Speedup only %.1fx (R: %.3fs, C++: %.3fs)",
                             speedup, time_r, time_cpp))
})

# ==============================================================================
# Test 10: Numerical Stability Tests
# ==============================================================================

test_that("MVNormal C++ maintains numerical stability", {
  skip_if_not(exists("mvnormal_likelihood_cpp"))

  # Test with very small eigenvalues (near-singular covariance)
  sigma_small_eig <- diag(c(1, 1e-8))
  x <- matrix(c(0, 0), ncol = 2)

  lik_small <- mvnormal_likelihood_cpp(x, c(0, 0), sigma_small_eig)
  expect_true(is.finite(lik_small))
  expect_true(lik_small > 0)

  # Test with very large values
  x_large <- matrix(c(1e6, 1e6), ncol = 2)
  lik_large <- mvnormal_likelihood_cpp(x_large, c(0, 0), diag(2))
  expect_true(is.finite(lik_large) || lik_large == 0)  # May underflow to 0

  # Test symmetry enforcement
  # Create slightly asymmetric matrix (due to numerical error)
  sigma_asym <- matrix(c(1, 0.5, 0.5 + 1e-15, 1), 2, 2)

  # Should still work (symmetry enforced internally)
  lik_asym <- mvnormal_likelihood_cpp(x, c(0, 0), sigma_asym)
  expect_true(is.finite(lik_asym))

  # Compare with symmetric version
  sigma_sym <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  lik_sym <- mvnormal_likelihood_cpp(x, c(0, 0), sigma_sym)
  expect_equal(lik_asym, lik_sym, tolerance = 1e-10)
})

# ==============================================================================
# Test 11: Consistency Tests Across Multiple Runs
# ==============================================================================

test_that("MVNormal C++ produces consistent results with same seed", {
  skip_if_not(can_use_cpp(DirichletProcessMvnormal(matrix(1), list())))

  data <- create_test_data(n = 50, d = 2, n_clusters = 2, seed = 999)
  priors <- create_test_priors(d = 2)

  # Run 1
  set_use_cpp(TRUE)
  set.seed(1234)
  dp1 <- DirichletProcessMvnormal(data, priors)
  dp1 <- Fit(dp1, 100, progressBar = FALSE)

  # Run 2 with same seed
  set.seed(1234)
  dp2 <- DirichletProcessMvnormal(data, priors)
  dp2 <- Fit(dp2, 100, progressBar = FALSE)

  # Results should be identical
  expect_equal(dp1$numberClusters, dp2$numberClusters)
  expect_equal(dp1$clusterLabels, dp2$clusterLabels)
  expect_equal(dp1$alphaChain, dp2$alphaChain)
  expect_equal(dp1$likelihoodChain, dp2$likelihoodChain)
})

# ==============================================================================
# Test 12: Integration with DP Methods
# ==============================================================================

test_that("MVNormal C++ integrates correctly with DP methods", {
  skip_if_not(can_use_cpp(DirichletProcessMvnormal(matrix(1), list())))

  data <- create_test_data(n = 100, d = 2, n_clusters = 2)
  priors <- create_test_priors(d = 2)

  set_use_cpp(TRUE)
  dp <- DirichletProcessMvnormal(data, priors)
  dp <- Fit(dp, 100, progressBar = FALSE)

  # Test posterior clusters
  post_clusters <- PosteriorClusters(dp)
  expect_type(post_clusters, "list")
  expect_true(length(post_clusters) > 0)

  # Test posterior function
  posterior_df <- PosteriorFrame(dp)
  expect_s3_class(posterior_df, "data.frame")
  expect_true(nrow(posterior_df) == nrow(data))
  expect_true("Mean" %in% names(posterior_df))
  expect_true("Variance" %in% names(posterior_df))

  # Test likelihood calculation
  lik <- Likelihood(dp)
  expect_type(lik, "double")
  expect_true(length(lik) == nrow(data))
  expect_true(all(lik > 0))

  # Test prediction for new data
  new_data <- create_test_data(n = 10, d = 2, n_clusters = 1)
  pred_probs <- PosteriorFunction(dp, new_data)
  expect_equal(dim(pred_probs), c(10, 100))  # 10 new points, 100 iterations
})

# Print summary
cat("\n=== MVNormal C++ Implementation Test Summary ===\n")
cat("This comprehensive test suite covers:\n")
cat("- Basic function availability\n")
cat("- Prior and posterior distributions\n")
cat("- Likelihood and predictive functions\n")
cat("- Cluster update algorithms\n")
cat("- Full MCMC integration\n")
cat("- Edge cases and numerical stability\n")
cat("- Performance benchmarks\n")
cat("- Consistency and reproducibility\n")
cat("- Integration with DP framework\n")
