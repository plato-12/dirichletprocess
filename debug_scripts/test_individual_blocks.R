# Test Individual Blocks from Comprehensive Test Suite
# Run each test_that block separately to isolate the crash

library(dirichletprocess)
library(testthat)
library(mvtnorm)

set_use_cpp(TRUE)

cat("=== Testing Individual Blocks ===\n")

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

# Test Block 1: MVNormal C++ Likelihood Function
cat("1. Testing Likelihood Function Block...\n")
tryCatch({
  test_that("MVNormal C++ Likelihood Function", {
    skip_if_not(cpp_available(), "C++ implementations not available")
    
    # Test data
    d <- 3
    n <- 50
    set.seed(42)
    
    # Create test parameters
    sigma <- diag(d) * 2
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
  cat("✓ PASSED: Likelihood Function Block\n")
}, error = function(e) {
  cat("✗ FAILED: Likelihood Function Block -", e$message, "\n")
})

# Test Block 2: MVNormal C++ Prior Draw Function  
cat("2. Testing Prior Draw Function Block...\n")
tryCatch({
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
      
      # Test prior draws
      n_draws <- 10
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
  cat("✓ PASSED: Prior Draw Function Block\n")
}, error = function(e) {
  cat("✗ FAILED: Prior Draw Function Block -", e$message, "\n")
})

# Test Block 3: MVNormal C++ Posterior Draw Function
cat("3. Testing Posterior Draw Function Block...\n")
tryCatch({
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
    
    # Test posterior draws
    n_draws <- 20
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
  cat("✓ PASSED: Posterior Draw Function Block\n")
}, error = function(e) {
  cat("✗ FAILED: Posterior Draw Function Block -", e$message, "\n")
})

# Test Block 4: MVNormal2 C++ Nonconjugate Implementation
cat("4. Testing MVNormal2 Nonconjugate Block...\n")
tryCatch({
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
  cat("✓ PASSED: MVNormal2 Nonconjugate Block\n")
}, error = function(e) {
  cat("✗ FAILED: MVNormal2 Nonconjugate Block -", e$message, "\n")
})

cat("=== Individual Block Testing Completed ===\n")
set_use_cpp(FALSE)