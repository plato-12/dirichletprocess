# Debug the benchmark framework error
library(dirichletprocess)

# Enable C++ implementations for better performance
set_use_cpp(TRUE)
enable_cpp_samplers()

cat("=== DEBUGGING BENCHMARK FRAMEWORK ERROR ===\n")

# Source the benchmark file to get the functions
source("benchmark/atime/benchmark-covariance-models-comprehensive.R")

# Test data generation
cat("Testing data generation...\n")
test_data_1d <- generate_benchmark_data(50, 1)
cat(sprintf("1D data dimensions: %s\n", paste(dim(test_data_1d), collapse = "x")))
cat(sprintf("1D data class: %s\n", class(test_data_1d)))

test_data_2d <- generate_benchmark_data(50, 2)
cat(sprintf("2D data dimensions: %s\n", paste(dim(test_data_2d), collapse = "x")))
cat(sprintf("2D data class: %s\n", class(test_data_2d)))

# Test parameter creation
cat("\nTesting parameter creation...\n")
tryCatch({
  prior_params_1d <- create_prior_parameters(1, "E")
  cat("1D parameters created successfully\n")
  print(str(prior_params_1d))
}, error = function(e) {
  cat(sprintf("1D parameter creation failed: %s\n", e$message))
})

tryCatch({
  prior_params_2d <- create_prior_parameters(2, "FULL")
  cat("2D parameters created successfully\n")
  print(str(prior_params_2d))
}, error = function(e) {
  cat(sprintf("2D parameter creation failed: %s\n", e$message))
})

# Test the collection function directly
cat("\nTesting performance metrics collection...\n")
tryCatch({
  # Test with simple 2D data
  test_data <- matrix(rnorm(100), ncol = 2)
  prior_params <- create_prior_parameters(2, "FULL")
  
  cat("Calling collect_performance_metrics...\n")
  result <- collect_performance_metrics("FULL", test_data, prior_params, 10, 5)
  
  cat("Performance metrics collected successfully\n")
  print(result)
  
}, error = function(e) {
  cat(sprintf("Performance metrics collection failed: %s\n", e$message))
  cat("Full error:\n")
  print(e)
})

cat("\n=== BENCHMARK DEBUG COMPLETE ===\n")