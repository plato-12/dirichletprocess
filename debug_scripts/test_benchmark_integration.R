# Test script for new benchmark integration functions

library(dirichletprocess)
devtools::load_all()

# Load atime if available
if (!require(atime, quietly = TRUE)) {
  cat("atime package not available\n")
  quit()
}

# Load mvtnorm for data generation
if (!require(mvtnorm, quietly = TRUE)) {
  cat("mvtnorm package not available\n")
  quit()
}

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test benchmark integration functions
test_benchmark_integration <- function() {
  cat("=== Testing Benchmark Integration Functions ===\n")
  
  # Source the new integration functions
  source("R/benchmark_integration.R")
  
  # Test 1: Test parameter preparation
  cat("\n1. Testing parameter preparation...\n")
  
  tryCatch({
    # Test valid cases
    params_2d_full <- prepare_benchmark_parameters(2, "FULL")
    params_1d_e <- prepare_benchmark_parameters(1, "E")
    params_3d_eii <- prepare_benchmark_parameters(3, "EII")
    
    cat("  Parameter preparation: SUCCESS\n")
    cat("    2D FULL:", length(params_2d_full$mu0), "dimensions\n")
    cat("    1D E model: SUCCESS\n")
    cat("    3D EII:", length(params_3d_eii$mu0), "dimensions\n")
    
  }, error = function(e) {
    cat("  Parameter preparation failed:", e$message, "\n")
  })
  
  # Test 2: Test data generation
  cat("\n2. Testing data generation...\n")
  
  tryCatch({
    # Test different data sizes and dimensions
    data_1d <- generate_benchmark_data(20, 1)
    data_2d <- generate_benchmark_data(30, 2)
    data_5d <- generate_benchmark_data(50, 5)
    
    cat("  Data generation: SUCCESS\n")
    cat("    1D data:", dim(data_1d), "\n")
    cat("    2D data:", dim(data_2d), "\n")
    cat("    5D data:", dim(data_5d), "\n")
    
  }, error = function(e) {
    cat("  Data generation failed:", e$message, "\n")
  })
  
  # Test 3: Test individual DP benchmark function
  cat("\n3. Testing individual DP benchmark...\n")
  
  tryCatch({
    # Test with different return formats
    test_data <- generate_benchmark_data(20, 2)
    
    result_metrics <- run_dp_benchmark(test_data, "FULL", mcmc_iterations = 3, return_format = "metrics")
    result_clusters <- run_dp_benchmark(test_data, "EII", mcmc_iterations = 3, return_format = "clusters")
    
    cat("  Individual DP benchmark: SUCCESS\n")
    cat("    Metrics format:", names(result_metrics), "\n")
    cat("    Clusters format:", result_clusters, "\n")
    
  }, error = function(e) {
    cat("  Individual DP benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 4: Test quick benchmark
  cat("\n4. Testing quick benchmark...\n")
  
  tryCatch({
    quick_result <- quick_benchmark_test()
    
    cat("  Quick benchmark: SUCCESS\n")
    cat("  Result structure:\n")
    str(quick_result, max.level = 2)
    
  }, error = function(e) {
    cat("  Quick benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 5: Test optimized atime benchmark
  cat("\n5. Testing optimized atime benchmark...\n")
  
  tryCatch({
    # Test with minimal settings
    optimized_result <- run_optimized_atime_benchmark(
      max_n = 30, 
      dimensions = 2, 
      models = c("FULL", "EII", "VII"), 
      mcmc_iter = 2, 
      repetitions = 2
    )
    
    cat("  Optimized atime benchmark: SUCCESS\n")
    cat("  Models tested:", unique(optimized_result$measurements$expr.name), "\n")
    cat("  Sample sizes tested:", unique(optimized_result$measurements$N), "\n")
    
  }, error = function(e) {
    cat("  Optimized atime benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 6: Test with univariate models
  cat("\n6. Testing univariate models...\n")
  
  tryCatch({
    # Test E and V models
    univariate_result <- run_optimized_atime_benchmark(
      max_n = 25, 
      dimensions = 1, 
      models = c("E", "V", "FULL"), 
      mcmc_iter = 2, 
      repetitions = 2
    )
    
    cat("  Univariate models benchmark: SUCCESS\n")
    cat("  Models tested:", unique(univariate_result$measurements$expr.name), "\n")
    
  }, error = function(e) {
    cat("  Univariate models benchmark failed:", e$message, "\n")
    print(e)
  })
  
  cat("\n=== Benchmark Integration Test Complete ===\n")
}

# Run the test
test_benchmark_integration()