# Minimal test to isolate C++ crash
library(dirichletprocess)

set_use_cpp(TRUE)

cat("Testing C++ availability...\n")
print(using_cpp())

# Test 1: Basic function existence
cat("Testing function existence...\n")
funcs_exist <- c(
  exists("mvnormal_likelihood_cpp"),
  exists("mvnormal_prior_draw_cpp"), 
  exists("mvnormal_posterior_draw_cpp")
)
print(funcs_exist)

if (!all(funcs_exist)) {
  cat("C++ functions not available\n")
  quit()
}

# Test 2: Very simple likelihood test
cat("Testing likelihood function...\n")
tryCatch({
  d <- 2
  mu <- c(0, 0)
  sigma <- diag(2)
  data <- matrix(c(1, 1), nrow = 1)
  
  lik <- mvnormal_likelihood_cpp(data, mu, sigma)
  cat("Likelihood test passed, result:", lik, "\n")
}, error = function(e) {
  cat("Likelihood test failed:", e$message, "\n")
})

# Test 3: Very simple prior draw
cat("Testing prior draw function...\n")
tryCatch({
  d <- 2
  priorParams <- list(
    mu0 = c(0, 0),
    kappa0 = 1,
    Lambda = diag(2),
    nu = 4
  )
  
  prior_draw <- mvnormal_prior_draw_cpp(priorParams, 1)
  cat("Prior draw test passed\n")
  print(str(prior_draw))
}, error = function(e) {
  cat("Prior draw test failed:", e$message, "\n")
})

# Test 4: Very simple posterior draw  
cat("Testing posterior draw function...\n")
tryCatch({
  d <- 2
  priorParams <- list(
    mu0 = c(0, 0),
    kappa0 = 1,
    Lambda = diag(2),
    nu = 4
  )
  
  test_data <- matrix(rnorm(4), 2, 2)
  
  posterior_draw <- mvnormal_posterior_draw_cpp(priorParams, test_data, 1)
  cat("Posterior draw test passed\n")
  print(str(posterior_draw))
}, error = function(e) {
  cat("Posterior draw test failed:", e$message, "\n")
})

cat("All tests completed without crash\n")