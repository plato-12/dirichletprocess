# Test the Weibull C++ implementation fix
library(dirichletprocess)

# Test function
test_weibull_cpp_fix <- function() {
  cat("Testing Weibull C++ implementation fix...\n")

  # Set seed for reproducibility
  set.seed(123)

  # Generate simple test data
  test_data <- rweibull(30, shape = 2, scale = 1)
  g0Priors <- c(10, 2, 4)  # phi, alpha0, beta0

  # Test 1: Basic functionality
  cat("\nTest 1: Basic C++ functionality\n")
  set_use_cpp(TRUE)

  tryCatch({
    dp_cpp <- DirichletProcessWeibull(test_data, g0Priors)
    dp_cpp <- Fit(dp_cpp, 10, progressBar = TRUE)
    cat("✓ Basic test passed - no crashes\n")
    cat("  Number of clusters:", dp_cpp$numberClusters, "\n")
    cat("  Alpha:", dp_cpp$alpha, "\n")
  }, error = function(e) {
    cat("✗ Basic test failed:", e$message, "\n")
    return(FALSE)
  })

  # Test 2: Longer run
  cat("\nTest 2: Extended run (50 iterations)\n")
  tryCatch({
    dp_cpp2 <- DirichletProcessWeibull(test_data, g0Priors)
    dp_cpp2 <- Fit(dp_cpp2, 50, progressBar = FALSE)
    cat("✓ Extended test passed\n")
    cat("  Final clusters:", dp_cpp2$numberClusters, "\n")
  }, error = function(e) {
    cat("✗ Extended test failed:", e$message, "\n")
    return(FALSE)
  })

  # Test 3: Edge cases
  cat("\nTest 3: Edge cases\n")

  # Small dataset
  small_data <- rweibull(5, shape = 1, scale = 1)
  tryCatch({
    dp_small <- DirichletProcessWeibull(small_data, g0Priors)
    dp_small <- Fit(dp_small, 20, progressBar = FALSE)
    cat("✓ Small dataset test passed\n")
  }, error = function(e) {
    cat("✗ Small dataset test failed:", e$message, "\n")
  })

  # Test 4: Compare with R implementation
  cat("\nTest 4: R vs C++ comparison\n")
  set.seed(456)

  set_use_cpp(FALSE)
  dp_r <- DirichletProcessWeibull(test_data, g0Priors)
  dp_r <- Fit(dp_r, 30, progressBar = FALSE)

  set_use_cpp(TRUE)
  set.seed(456)
  dp_cpp3 <- DirichletProcessWeibull(test_data, g0Priors)
  dp_cpp3 <- Fit(dp_cpp3, 30, progressBar = FALSE)

  cat("  R implementation - Clusters:", dp_r$numberClusters, "\n")
  cat("  C++ implementation - Clusters:", dp_cpp3$numberClusters, "\n")
  cat("  Difference:", abs(dp_r$numberClusters - dp_cpp3$numberClusters), "\n")

  if (abs(dp_r$numberClusters - dp_cpp3$numberClusters) <= 2) {
    cat("✓ Results are reasonably similar\n")
  } else {
    cat("⚠ Results differ significantly\n")
  }

  cat("\nAll tests completed!\n")
  return(TRUE)
}

# Run the test
test_weibull_cpp_fix()
