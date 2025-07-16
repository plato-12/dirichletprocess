# Test the optimized benchmark system

library(dirichletprocess)
devtools::load_all()

# Load required packages
if (!require(atime, quietly = TRUE)) {
  cat("atime package not available\n")
  quit()
}
if (!require(mvtnorm, quietly = TRUE)) {
  cat("mvtnorm package not available\n")
  quit()
}

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test optimized benchmark
test_optimized_benchmark <- function() {
  cat("=== Testing Optimized Benchmark System ===\n")
  
  # Source the optimized benchmark
  source("benchmark/atime/benchmark-covariance-models-optimized.R")
  
  # Test 1: Quick test of all models
  cat("\n1. Running quick test of all models...\n")
  
  tryCatch({
    quick_results <- quick_test_all_models()
    
    cat("  Quick test: SUCCESS\n")
    cat("  Univariate models tested:", nrow(quick_results$univariate$measurements), "measurements\n")
    cat("  Multivariate models tested:", nrow(quick_results$multivariate$measurements), "measurements\n")
    
  }, error = function(e) {
    cat("  Quick test failed:", e$message, "\n")
    print(e)
  })
  
  # Test 2: Standard benchmark
  cat("\n2. Running standard benchmark...\n")
  
  tryCatch({
    standard_result <- run_fast_covariance_benchmark(STANDARD_CONFIG, dimensions = 2)
    analysis <- print_benchmark_summary(standard_result)
    
    cat("  Standard benchmark: SUCCESS\n")
    
  }, error = function(e) {
    cat("  Standard benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 3: Univariate-only benchmark
  cat("\n3. Running univariate-only benchmark...\n")
  
  tryCatch({
    univariate_result <- run_univariate_benchmark(FAST_CONFIG)
    
    cat("  Univariate benchmark: SUCCESS\n")
    cat("  Models tested:", paste(unique(univariate_result$measurements$expr.name), collapse = ", "), "\n")
    
  }, error = function(e) {
    cat("  Univariate benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 4: Multivariate-only benchmark
  cat("\n4. Running multivariate-only benchmark...\n")
  
  tryCatch({
    multivariate_result <- run_multivariate_benchmark(FAST_CONFIG, dimensions = 3)
    
    cat("  Multivariate benchmark: SUCCESS\n")
    cat("  Models tested:", paste(unique(multivariate_result$measurements$expr.name), collapse = ", "), "\n")
    
  }, error = function(e) {
    cat("  Multivariate benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 5: Test comprehensive multi-dimensional
  cat("\n5. Running comprehensive multi-dimensional benchmark...\n")
  
  tryCatch({
    comprehensive_results <- run_comprehensive_benchmark(FAST_CONFIG, c(1, 2))
    
    cat("  Comprehensive benchmark: SUCCESS\n")
    cat("  Dimensions tested:", names(comprehensive_results), "\n")
    
  }, error = function(e) {
    cat("  Comprehensive benchmark failed:", e$message, "\n")
    print(e)
  })
  
  cat("\n=== Optimized Benchmark Test Complete ===\n")
}

# Run the test
test_optimized_benchmark()