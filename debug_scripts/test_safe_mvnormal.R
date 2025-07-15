# Safe testing version of MVNormal C++ comprehensive tests
# This version adds error handling to isolate crash points

library(dirichletprocess)
library(testthat)
library(mvtnorm)

set_use_cpp(TRUE)

cat("=== Starting Safe MVNormal C++ Tests ===\n")

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

# Test 1: Basic C++ availability
safe_test("C++ Availability Check", function() {
  stopifnot(using_cpp())
  stopifnot(exists("mvnormal_likelihood_cpp"))
  stopifnot(exists("mvnormal_prior_draw_cpp"))
  stopifnot(exists("mvnormal_posterior_draw_cpp"))
})

# Test 2: Simple likelihood 
safe_test("Simple Likelihood Test", function() {
  d <- 2
  mu <- c(0, 0)
  sigma <- diag(2)
  data <- matrix(c(1, 1), nrow = 1)
  
  lik <- mvnormal_likelihood_cpp(data, mu, sigma)
  stopifnot(is.finite(lik))
  stopifnot(lik > 0)
})

# Test 3: Prior draw with minimal parameters
safe_test("Minimal Prior Draw Test", function() {
  d <- 2
  priorParams <- list(
    mu0 = c(0, 0),
    kappa0 = 1,
    Lambda = diag(2),
    nu = 4
  )
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 1)
  stopifnot(is.list(prior_draw))
  stopifnot(all(c("mu", "sig") %in% names(prior_draw)))
  
  # Check dimensions
  cat("Prior draw mu dimensions:", dim(prior_draw$mu), "\n")
  cat("Prior draw sig dimensions:", dim(prior_draw$sig), "\n")
  
  stopifnot(all(is.finite(prior_draw$mu)))
  stopifnot(all(is.finite(prior_draw$sig)))
})

# Test 4: Prior draw with multiple draws
safe_test("Multiple Prior Draws Test", function() {
  d <- 2
  priorParams <- list(
    mu0 = c(0, 0),
    kappa0 = 1,
    Lambda = diag(2),
    nu = 4
  )
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 3)
  stopifnot(is.list(prior_draw))
  
  cat("Multiple draws mu dimensions:", dim(prior_draw$mu), "\n")
  cat("Multiple draws sig dimensions:", dim(prior_draw$sig), "\n")
  
  stopifnot(all(is.finite(prior_draw$mu)))
  stopifnot(all(is.finite(prior_draw$sig)))
})

# Test 5: Higher dimensional prior draw
safe_test("Higher Dimensional Prior Draw", function() {
  d <- 3
  priorParams <- list(
    mu0 = rep(0, d),
    kappa0 = 1,
    Lambda = diag(d),
    nu = d + 2
  )
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 2)
  
  cat("Higher dim mu dimensions:", dim(prior_draw$mu), "\n")
  cat("Higher dim sig dimensions:", dim(prior_draw$sig), "\n")
  
  stopifnot(all(is.finite(prior_draw$mu)))
  stopifnot(all(is.finite(prior_draw$sig)))
})

# Test 6: Posterior draw
safe_test("Posterior Draw Test", function() {
  d <- 2
  priorParams <- list(
    mu0 = c(0, 0),
    kappa0 = 1,
    Lambda = diag(2),
    nu = 4
  )
  
  test_data <- matrix(rnorm(6), 3, 2)
  
  posterior_draw <- mvnormal_posterior_draw_cpp(priorParams, test_data, 2)
  
  cat("Posterior mu dimensions:", dim(posterior_draw$mu), "\n")
  cat("Posterior sig dimensions:", dim(posterior_draw$sig), "\n")
  
  stopifnot(all(is.finite(posterior_draw$mu)))
  stopifnot(all(is.finite(posterior_draw$sig)))
})

# Test 7: MVNormal2 functions
safe_test("MVNormal2 Functions Test", function() {
  if (!exists("mvnormal2_prior_draw_cpp")) {
    cat("MVNormal2 functions not available, skipping\n")
    return(TRUE)
  }
  
  d <- 2
  priorParams <- list(
    mu0 = matrix(rep(0, d), nrow = 1),
    sigma0 = diag(d),
    phi0 = diag(d) * 2,
    nu0 = d + 2
  )
  
  prior_draw <- mvnormal2_prior_draw_cpp(priorParams, 2)
  
  cat("MVNormal2 mu dimensions:", dim(prior_draw$mu), "\n")
  cat("MVNormal2 sig dimensions:", dim(prior_draw$sig), "\n")
  
  stopifnot(all(is.finite(prior_draw$mu)))
  stopifnot(all(is.finite(prior_draw$sig)))
})

# Test 8: DirichletProcessMvnormal2 creation
safe_test("DirichletProcessMvnormal2 Creation", function() {
  set.seed(123)
  test_data <- matrix(rnorm(20), 10, 2)
  
  dp <- DirichletProcessMvnormal2(test_data)
  
  stopifnot(inherits(dp, "nonconjugate"))
  stopifnot(inherits(dp$mixingDistribution, "mvnormal2"))
})

cat("=== Safe MVNormal C++ Tests Completed ===\n")
set_use_cpp(FALSE)