# Test testthat Integration Issues
# This isolates testthat-specific problems

library(dirichletprocess)
library(testthat)
library(mvtnorm)

set_use_cpp(TRUE)

cat("=== Testing testthat Integration ===\n")

# Test individual expect_* functions to find the problematic one

cat("1. Testing basic expect_true...\n")
tryCatch({
  expect_true(TRUE)
  cat("✓ expect_true works\n")
}, error = function(e) cat("✗ expect_true failed:", e$message, "\n"))

cat("2. Testing expect_equal with simple values...\n")
tryCatch({
  expect_equal(1, 1)
  cat("✓ expect_equal works\n")
}, error = function(e) cat("✗ expect_equal failed:", e$message, "\n"))

cat("3. Testing C++ function call in testthat context...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
  cat("✓ C++ function call works\n")
}, error = function(e) cat("✗ C++ function call failed:", e$message, "\n"))

cat("4. Testing expect_true with C++ result...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
  expect_true(is.list(prior_draw))
  cat("✓ expect_true with C++ result works\n")
}, error = function(e) cat("✗ expect_true with C++ result failed:", e$message, "\n"))

cat("5. Testing expect_equal with dimensions...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  n_draws <- 5
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, n_draws)
  
  cat("  Checking dim result:", dim(prior_draw$mu), "\n")
  cat("  Expected:", c(1, d, n_draws), "\n")
  
  expect_equal(dim(prior_draw$mu), c(1, d, n_draws))
  cat("✓ expect_equal with dimensions works\n")
}, error = function(e) cat("✗ expect_equal with dimensions failed:", e$message, "\n"))

cat("6. Testing expect_true with all() and is.finite()...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
  
  cat("  Checking is.finite on mu...\n")
  finite_check <- all(is.finite(prior_draw$mu))
  cat("  is.finite result:", finite_check, "\n")
  
  expect_true(finite_check)
  cat("✓ expect_true with all(is.finite()) works\n")
}, error = function(e) cat("✗ expect_true with all(is.finite()) failed:", e$message, "\n"))

cat("7. Testing eigenvalue computation in testthat...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
  
  for (i in 1:5) {
    cat("  Computing eigenvalues for draw", i, "\n")
    eigenvals <- eigen(prior_draw$sig[, , i])$values
    cat("  Eigenvalues:", eigenvals, "\n")
    
    expect_true(all(eigenvals > 1e-10))
  }
  cat("✓ Eigenvalue computation in testthat works\n")
}, error = function(e) cat("✗ Eigenvalue computation in testthat failed:", e$message, "\n"))

cat("8. Testing apply() operation in testthat...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 10)
  
  cat("  Array structure before apply:\n")
  str(prior_draw$mu)
  
  cat("  Performing apply operation...\n")
  mu_means <- apply(prior_draw$mu[1, , ], 2, mean)
  cat("  Apply result:", head(mu_means), "...\n")
  
  expect_equal(length(mu_means), d)
  cat("✓ apply() operation in testthat works\n")
}, error = function(e) cat("✗ apply() operation in testthat failed:", e$message, "\n"))

cat("9. Testing with tolerance in expect_equal...\n")
tryCatch({
  d <- 2
  mu0 <- rep(0, d)
  kappa0 <- 1
  Lambda <- diag(d)
  nu <- d + 2
  priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 10)
  mu_means <- apply(prior_draw$mu[1, , ], 2, mean)
  
  expect_equal(mu_means, mu0, tolerance = 0.5)
  cat("✓ expect_equal with tolerance works\n")
}, error = function(e) cat("✗ expect_equal with tolerance failed:", e$message, "\n"))

cat("10. Testing full test_that wrapper...\n")
tryCatch({
  test_that("Test wrapper works", {
    d <- 2
    mu0 <- rep(0, d)
    kappa0 <- 1
    Lambda <- diag(d)
    nu <- d + 2
    priorParams <- list(mu0 = mu0, kappa0 = kappa0, Lambda = Lambda, nu = nu)
    
    prior_draw <- mvnormal_prior_draw_cpp(priorParams, 5)
    
    expect_true(is.list(prior_draw))
    expect_equal(dim(prior_draw$mu), c(1, d, 5))
    expect_true(all(is.finite(prior_draw$mu)))
  })
  cat("✓ Full test_that wrapper works\n")
}, error = function(e) cat("✗ Full test_that wrapper failed:", e$message, "\n"))

cat("=== testthat Integration Test Completed ===\n")
set_use_cpp(FALSE)