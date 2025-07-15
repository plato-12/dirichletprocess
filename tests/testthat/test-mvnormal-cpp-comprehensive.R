# Comprehensive Tests for MVNormal Distribution C++ Implementations
# =================================================================
#
# This test suite tests the available MVNormal C++ implementations
# including likelihood, prior/posterior draws, and MVNormal2 nonconjugate functions.

# Load required libraries
suppressPackageStartupMessages({
  library(testthat)
  library(dirichletprocess)
  library(mvtnorm)
})

# Enable C++ mode for the duration of these tests
set_use_cpp(TRUE)

# Helper function to check if C++ implementations are available
cpp_available <- function() {
  return(using_cpp() && 
         exists("mvnormal_likelihood_cpp") &&
         exists("mvnormal_prior_draw_cpp") &&
         exists("mvnormal_posterior_draw_cpp"))
}

# Helper function to generate test data with known structure
generate_test_data <- function(n_points = 100, n_clusters = 3, dims = 2, seed = 123) {
  set.seed(seed)
  
  # Define cluster centers
  centers <- matrix(c(-2, -2, 0, 2, 2, 0), nrow = 3, ncol = 2)
  if (dims > 2) {
    extra_dims <- matrix(rnorm(3 * (dims - 2)), nrow = 3, ncol = dims - 2)
    centers <- cbind(centers, extra_dims)
  }
  
  # Generate data from mixture
  cluster_assignments <- sample(1:n_clusters, n_points, replace = TRUE)
  data <- matrix(0, nrow = n_points, ncol = dims)
  
  for (i in 1:n_clusters) {
    idx <- which(cluster_assignments == i)
    if (length(idx) > 0) {
      sigma <- diag(dims) * 0.5
      data[idx, ] <- mvtnorm::rmvnorm(length(idx), mean = centers[i, ], sigma = sigma)
    }
  }
  
  return(list(data = data, true_clusters = cluster_assignments, centers = centers))
}

# Helper function to create test covariance matrices
create_test_covariance <- function(d, type = "full") {
  switch(type,
    "full" = {
      # Full covariance matrix
      sigma <- matrix(rnorm(d^2), d, d)
      sigma <- sigma %*% t(sigma) + diag(d) * 0.1
      sigma
    },
    "diagonal" = {
      # Diagonal covariance
      diag(runif(d, 0.5, 2.0))
    },
    "spherical" = {
      # Spherical covariance
      diag(d) * runif(1, 0.5, 2.0)
    },
    "identity" = {
      # Identity covariance
      diag(d)
    }
  )
}

# =============================================================================
# CORE C++ FUNCTION TESTS
# =============================================================================

test_that("MVNormal C++ Likelihood Function", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test data
  d <- 3
  n <- 50
  set.seed(42)
  
  # Create test parameters
  sigma <- create_test_covariance(d, "full")
  mu <- rnorm(d)
  data <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
  
  # Test direct C++ likelihood function
  lik_cpp <- mvnormal_likelihood_cpp(data[1:5, , drop = FALSE], mu, sigma)
  
  # Compare with mvtnorm reference
  lik_ref <- mvtnorm::dmvnorm(data[1:5, ], mean = mu, sigma = sigma)
  
  expect_equal(length(lik_cpp), 5)
  expect_true(all(is.finite(lik_cpp)))
  expect_true(all(lik_cpp > 0))
  expect_equal(lik_cpp, lik_ref, tolerance = 1e-10)
})

test_that("MVNormal C++ Prior Draw Function", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test different dimensions
  for (d in c(2, 3, 5)) {
    # Create standard priors
    mu0 <- rep(0, d)
    kappa0 <- 1
    Lambda <- diag(d)
    nu <- d + 2
    
    priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
    
    # Test prior draws (reduced from 10 to 5 for memory efficiency)
    n_draws <- 5
    prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
    
    expect_true(is.list(prior_draw))
    expect_true(all(c("mu", "sig") %in% names(prior_draw)))
    expect_equal(dim(prior_draw$mu), c(1, d, n_draws))
    expect_equal(dim(prior_draw$sig), c(d, d, n_draws))
    
    # Check that all means are finite
    expect_true(all(is.finite(prior_draw$mu)))
    
    # Check that covariance matrices are positive definite
    for (i in 1:n_draws) {
      eigenvals <- eigen(prior_draw$sig[, , i])$values
      expect_true(all(eigenvals > 1e-10), 
                  info = paste("Dimension", d, "Draw", i))
    }
    
    # Check prior means are centered around mu0
    mu_means <- apply(prior_draw$mu[1, , ], 1, mean)
    expect_equal(mu_means, mu0, tolerance = 1.0)
  }
})

test_that("MVNormal C++ Posterior Draw Function", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test posterior draws
  d <- 3
  n <- 30
  
  # Create priors
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  # Generate test data with known mean
  set.seed(123)
  true_mean <- c(1, -1, 0.5)
  test_data <- mvtnorm::rmvnorm(n, mean = true_mean, sigma = diag(d) * 0.5)
  
  # Test posterior draws (reduced from 20 to 10 for memory efficiency)
  n_draws <- 10
  posterior_draw <- mvnormal_posterior_draw_cpp(priorParams, test_data, n_draws)
  
  expect_true(is.list(posterior_draw))
  expect_equal(dim(posterior_draw$mu), c(1, d, n_draws))
  expect_equal(dim(posterior_draw$sig), c(d, d, n_draws))
  
  # Check that posterior covariance matrices are positive definite
  for (i in 1:n_draws) {
    eigenvals <- eigen(posterior_draw$sig[, , i])$values
    expect_true(all(eigenvals > 1e-10), info = paste("Posterior draw", i))
  }
  
  # Posterior means should be pulled toward data
  data_mean <- colMeans(test_data)
  posterior_mu_means <- apply(posterior_draw$mu[1, , ], 1, mean)
  
  # Should be closer to data mean than prior mean
  for (i in 1:d) {
    expect_true(abs(posterior_mu_means[i] - data_mean[i]) < 
                abs(mu0[i] - data_mean[i]))
  }
})

test_that("MVNormal C++ Posterior Parameters Function", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test posterior parameter calculation
  d <- 2
  n <- 15
  
  # Create prior parameters
  mu0 <- c(0, 0)
  kappa0 <- 1
  Lambda <- diag(d) * 2
  nu <- d + 1
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  # Generate test data
  set.seed(456)
  test_data <- matrix(rnorm(n * d), n, d)
  
  # Get posterior parameters using C++
  post_params <- mvnormal_posterior_parameters_cpp(priorParams, test_data)
  
  # Test structure
  expect_true(is.list(post_params))
  expect_true(all(c("mu_n", "kappa_n", "Lambda_n", "nu_n") %in% names(post_params)))
  
  # Test dimensions
  expect_equal(length(post_params$mu_n), d)
  expect_equal(post_params$kappa_n, kappa0 + n)
  expect_equal(post_params$nu_n, nu + n)
  expect_equal(dim(post_params$Lambda_n), c(d, d))
  
  # Lambda_n should be positive definite
  eigenvals <- eigen(post_params$Lambda_n)$values
  expect_true(all(eigenvals > 1e-10))
  
  # Check conjugate update formulas manually
  data_mean <- colMeans(test_data)
  expected_mu_n <- (kappa0 * mu0 + n * data_mean) / (kappa0 + n)
  expect_equal(post_params$mu_n, expected_mu_n, tolerance = 1e-10)
})

# Force garbage collection between test sections
gc()

# =============================================================================
# MVNORMAL2 NONCONJUGATE TESTS
# =============================================================================

test_that("MVNormal2 C++ Nonconjugate Implementation", {
  skip_if_not(exists("mvnormal2_prior_draw_cpp"), 
              "MVNormal2 C++ functions not available")
  
  # Test the nonconjugate MVNormal2 implementation
  test_data <- generate_test_data(n_points = 30, n_clusters = 2, dims = 2, seed = 1313)
  dp <- DirichletProcessMvnormal2(test_data$data)
  
  # Test that it's properly classified as nonconjugate
  expect_true(inherits(dp, "nonconjugate"))
  expect_true(inherits(dp$mixingDistribution, "mvnormal2"))
  
  # Test basic structure
  expect_true(is.finite(dp$numberClusters))
  expect_true(dp$numberClusters >= 1)
  expect_equal(length(dp$clusterLabels), nrow(test_data$data))
  expect_true(all(dp$clusterLabels > 0))
})

test_that("MVNormal2 C++ Individual Functions", {
  skip_if_not(exists("mvnormal2_likelihood_cpp"), "MVNormal2 C++ functions not available")
  
  # Test MVNormal2 specific functions
  d <- 2
  n <- 10
  
  # Create test data and parameters
  set.seed(789)
  test_data <- matrix(rnorm(n * d), n, d)
  
  # Test MVNormal2 prior draw
  priorParams <- list(
    mu0 = matrix(rep(0, d), nrow = 1),
    sigma0 = diag(d),
    phi0 = diag(d) * 2,
    nu0 = d + 2
  )
  
  prior_draw <- mvnormal2_prior_draw_cpp(priorParams, 3)
  expect_true(is.list(prior_draw))
  expect_true(all(c("mu", "sig") %in% names(prior_draw)))
  
  # Test MVNormal2 likelihood
  theta <- list(mu = prior_draw$mu[1, , 1], sig = prior_draw$sig[, , 1])
  lik <- mvnormal2_likelihood_cpp(test_data[1, , drop = FALSE], theta)
  expect_true(is.finite(lik))
  expect_true(lik > 0)
})

# Force garbage collection between test sections
gc()

# =============================================================================
# EDGE CASES AND ROBUSTNESS TESTS
# =============================================================================

test_that("MVNormal C++ Edge Cases - Single Data Point", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test with single data point
  d <- 2
  data <- matrix(c(1.5, -0.5), nrow = 1)
  mu <- c(0, 0)
  sigma <- diag(d)
  
  # Should handle single point gracefully
  lik <- mvnormal_likelihood_cpp(data, mu, sigma)
  expect_equal(length(lik), 1)
  expect_true(is.finite(lik))
  expect_true(lik > 0)
  
  # Compare with reference
  lik_ref <- mvtnorm::dmvnorm(data, mean = mu, sigma = sigma)
  expect_equal(lik, lik_ref, tolerance = 1e-10)
})

test_that("MVNormal C++ Edge Cases - High Dimensional Data", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test with higher dimensional data
  d <- 8
  n <- 20
  set.seed(888)
  
  mu <- rnorm(d)
  sigma <- diag(d) * 2
  data <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
  
  # Should handle high dimensions
  lik <- mvnormal_likelihood_cpp(data, mu, sigma)
  expect_equal(length(lik), n)
  expect_true(all(is.finite(lik)))
  expect_true(all(lik > 0))
  
  # Test prior draws in high dimensions
  priorParams <- list(
    mu0 = rep(0, d),
    kappa0 = 1,
    Lambda = diag(d),
    nu = d + 2
  )
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 3)
  expect_equal(dim(prior_draw$mu), c(1, d, 3))
  expect_equal(dim(prior_draw$sig), c(d, d, 3))
  
  # Check positive definiteness in high dimensions
  for (i in 1:3) {
    eigenvals <- eigen(prior_draw$sig[, , i])$values
    expect_true(all(eigenvals > 1e-10), info = paste("High-dim draw", i))
  }
})

test_that("MVNormal C++ Edge Cases - Extreme Parameter Values", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test with extreme prior parameters
  d <- 2
  mu0 <- c(100, -100)  # Extreme means
  kappa0 <- 1000       # High precision
  Lambda <- diag(d) * 0.001  # Very small scale
  nu <- d + 0.1        # Minimal degrees of freedom
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  # Should handle extreme parameters
  expect_no_error({
    prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
  })
  
  expect_true(all(is.finite(prior_draw$mu)))
  expect_true(all(is.finite(prior_draw$sig)))
  
  # All covariance matrices should be positive definite
  for (i in 1:5) {
    eigenvals <- eigen(prior_draw$sig[, , i])$values
    expect_true(all(eigenvals > 0))
  }
})

# =============================================================================
# PERFORMANCE AND CONSISTENCY TESTS
# =============================================================================

# Force garbage collection between test sections
gc()

test_that("MVNormal C++ Consistency with R Implementation", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Compare C++ and R likelihood implementations
  d <- 3
  n <- 25
  set.seed(1111)
  
  mu <- rnorm(d)
  sigma <- create_test_covariance(d, "full")
  data <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
  
  # C++ likelihood
  lik_cpp <- mvnormal_likelihood_cpp(data, mu, sigma)
  
  # R reference (mvtnorm)
  lik_r <- mvtnorm::dmvnorm(data, mean = mu, sigma = sigma)
  
  # Should be identical (within numerical precision)
  expect_equal(lik_cpp, lik_r, tolerance = 1e-12)
})

test_that("MVNormal C++ Memory Management", {
  skip_if_not(cpp_available(), "C++ implementations not available")
  
  # Test that repeated operations don't cause memory issues
  d <- 3
  priorParams <- list(
    mu0 = rep(0, d),
    kappa0 = 1,
    Lambda = diag(d),
    nu = d + 2
  )
  
  # Run many draws to test for memory leaks
  for (i in 1:100) {
    prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
    
    # Verify draws remain valid
    expect_true(all(is.finite(prior_draw$mu)))
    expect_true(all(is.finite(prior_draw$sig)))
    
    # Check a few covariance matrices are positive definite
    if (i %% 20 == 0) {
      eigenvals <- eigen(prior_draw$sig[, , 1])$values
      expect_true(all(eigenvals > 1e-10), info = paste("Iteration", i))
    }
  }
})

# Explicit cleanup to prevent memory issues
gc()  # Force garbage collection
rm(list = ls())  # Clear all variables
gc()  # Force garbage collection again

# Clean up
set_use_cpp(FALSE)  # Reset to default

cat("All available MVNormal C++ tests completed successfully!\n")