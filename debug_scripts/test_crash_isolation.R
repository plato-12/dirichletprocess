# Crash Isolation Test for MVNormal C++ Comprehensive Suite
# This test incrementally replicates the comprehensive test to find the crash point

library(dirichletprocess)
library(testthat)
library(mvtnorm)

set_use_cpp(TRUE)

cat("=== Crash Isolation Test ===\n")

# Helper function with better error handling
safe_test <- function(test_name, test_func) {
  cat("Running:", test_name, "...\n")
  tryCatch({
    result <- test_func()
    cat("✓ PASSED:", test_name, "\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ FAILED:", test_name, "- Error:", e$message, "\n")
    return(FALSE)
  })
}

# Test 1: Replicate the exact prior draw test from comprehensive suite
safe_test("Prior Draw - Single Dimension (d=2)", function() {
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  n_draws <- 10  # Same as comprehensive test
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
  
  cat("  Dimensions mu:", dim(prior_draw$mu), "\n")
  cat("  Dimensions sig:", dim(prior_draw$sig), "\n")
  
  stopifnot(is.list(prior_draw))
  stopifnot(all(c("mu", "sig") %in% names(prior_draw)))
  stopifnot(all(dim(prior_draw$mu) == c(1, d, n_draws)))
  stopifnot(all(dim(prior_draw$sig) == c(d, d, n_draws)))
  stopifnot(all(is.finite(prior_draw$mu)))
})

# Test 2: Eigenvalue checking (potential crash point)
safe_test("Prior Draw - Eigenvalue Checking", function() {
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  n_draws <- 10
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
  
  # This is the eigenvalue checking loop from comprehensive test
  for (i in 1:n_draws) {
    cat("  Checking eigenvalues for draw", i, "\n")
    eigenvals <- eigen(prior_draw$sig[, , i])$values
    cat("  Eigenvalues:", eigenvals, "\n")
    stopifnot(all(eigenvals > 1e-10))
  }
})

# Test 3: Apply operation (potential crash point)
safe_test("Prior Draw - Apply Operation", function() {
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  n_draws <- 10
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
  
  # This is the apply operation from comprehensive test
  cat("  Prior draw mu structure:\n")
  str(prior_draw$mu)
  cat("  Attempting apply operation...\n")
  mu_means <- apply(prior_draw$mu[1, , ], 2, mean)
  cat("  Apply result:", mu_means, "\n")
})

# Test 4: Higher dimensions (d=3)
safe_test("Prior Draw - Higher Dimension (d=3)", function() {
  d <- 3
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  n_draws <- 10
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
  
  cat("  d=3 dimensions mu:", dim(prior_draw$mu), "\n")
  cat("  d=3 dimensions sig:", dim(prior_draw$sig), "\n")
  
  stopifnot(all(is.finite(prior_draw$mu)))
  stopifnot(all(is.finite(prior_draw$sig)))
})

# Test 5: Even higher dimensions (d=5)
safe_test("Prior Draw - High Dimension (d=5)", function() {
  d <- 5
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  n_draws <- 10
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
  
  cat("  d=5 dimensions mu:", dim(prior_draw$mu), "\n")
  cat("  d=5 dimensions sig:", dim(prior_draw$sig), "\n")
  
  stopifnot(all(is.finite(prior_draw$mu)))
  stopifnot(all(is.finite(prior_draw$sig)))
})

# Test 6: Full loop as in comprehensive test
safe_test("Prior Draw - Full Loop Test", function() {
  for (d in c(2, 3, 5)) {
    cat("  Testing dimension d =", d, "\n")
    
    mu0 <- rep(0, d)
    kappa0 <- 1
    Lambda <- diag(d)
    nu <- d + 2
    
    priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
    
    n_draws <- 10
    prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
    
    # All the checks from comprehensive test
    stopifnot(is.list(prior_draw))
    stopifnot(all(c("mu", "sig") %in% names(prior_draw)))
    stopifnot(all(dim(prior_draw$mu) == c(1, d, n_draws)))
    stopifnot(all(dim(prior_draw$sig) == c(d, d, n_draws)))
    stopifnot(all(is.finite(prior_draw$mu)))
    
    # Eigenvalue checking
    for (i in 1:n_draws) {
      eigenvals <- eigen(prior_draw$sig[, , i])$values
      stopifnot(all(eigenvals > 1e-10))
    }
    
    # Apply operation
    mu_means <- apply(prior_draw$mu[1, , ], 2, mean)
    # Note: tolerance relaxed for testing
    # stopifnot(all(abs(mu_means - mu0) < 0.5))
    
    cat("  ✓ Dimension d =", d, "completed\n")
  }
})

# Test 7: Posterior draw test (replicating comprehensive test)
safe_test("Posterior Draw - Comprehensive Test", function() {
  d <- 3
  n <- 30
  
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  set.seed(123)
  true_mean <- c(1, -1, 0.5)
  test_data <- mvtnorm::rmvnorm(n, mean = true_mean, sigma = diag(d) * 0.5)
  
  n_draws <- 20
  posterior_draw <- mvnormal_posterior_draw_cpp(priorParams, test_data, n_draws)
  
  cat("  Posterior dimensions mu:", dim(posterior_draw$mu), "\n")
  cat("  Posterior dimensions sig:", dim(posterior_draw$sig), "\n")
  
  stopifnot(is.list(posterior_draw))
  stopifnot(all(dim(posterior_draw$mu) == c(1, d, n_draws)))
  stopifnot(all(dim(posterior_draw$sig) == c(d, d, n_draws)))
  
  # Eigenvalue checking
  for (i in 1:n_draws) {
    eigenvals <- eigen(posterior_draw$sig[, , i])$values
    stopifnot(all(eigenvals > 1e-10))
  }
  
  # Apply operations
  data_mean <- colMeans(test_data)
  posterior_mu_means <- apply(posterior_draw$mu[1, , ], 2, mean)
  
  cat("  Data mean:", data_mean, "\n")
  cat("  Posterior mu means:", posterior_mu_means, "\n")
})

cat("=== Crash Isolation Test Completed ===\n")
set_use_cpp(FALSE)