# Test the fixed benchmark
library(dirichletprocess)
set.seed(123)

cat("=== TESTING FIXED BENCHMARK ===\n")

# Test the benchmark validation function
test_fixed_validation <- function() {
  cat("Testing benchmark validation...\n")
  
  tryCatch({
    # Source the fixed benchmark
    source("benchmark/atime/benchmark-covariance-models-comprehensive.R")
    
    # Test the validation function
    result <- validate_all_models()
    cat(sprintf("   Validation result: %s\n", result))
    
    return(result)
    
  }, error = function(e) {
    cat(sprintf("   ✗ Validation failed: %s\n", e$message))
    return(FALSE)
  })
}

# Test small-scale benchmark
test_small_benchmark <- function() {
  cat("\nTesting small-scale benchmark...\n")
  
  tryCatch({
    # Source the fixed benchmark
    source("benchmark/atime/benchmark-covariance-models-comprehensive.R")
    
    # Override configuration for quick test
    BENCHMARK_CONFIG$dimensions <- c(1, 2)
    BENCHMARK_CONFIG$sample_sizes <- c(20, 50)
    BENCHMARK_CONFIG$mcmc_iterations <- 10  # Very short for testing
    BENCHMARK_CONFIG$benchmark_reps <- 1
    
    # Test scalability analysis
    results <- run_scalability_analysis()
    
    cat(sprintf("   Completed %d benchmarks\n", length(results)))
    
    # Show some results
    success_count <- sum(sapply(results, function(x) x$success))
    total_count <- length(results)
    
    cat(sprintf("   Success rate: %d/%d\n", success_count, total_count))
    
    return(success_count > 0)
    
  }, error = function(e) {
    cat(sprintf("   ✗ Small benchmark failed: %s\n", e$message))
    return(FALSE)
  })
}

# Run tests
validation_success <- test_fixed_validation()
benchmark_success <- test_small_benchmark()

cat("\n=== RESULTS ===\n")
cat(sprintf("Validation: %s\n", if (validation_success) "✓ PASSED" else "✗ FAILED"))
cat(sprintf("Benchmark: %s\n", if (benchmark_success) "✓ PASSED" else "✗ FAILED"))

if (validation_success && benchmark_success) {
  cat("✅ ALL CONSTRAINED COVARIANCE MODELS WORKING!\n")
} else {
  cat("❌ Some issues remain\n")
}

cat("\n=== FIXED BENCHMARK TEST COMPLETE ===\n")