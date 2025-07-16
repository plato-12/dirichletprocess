# Test benchmark framework fixes
library(dirichletprocess)
set_use_cpp(TRUE)

# Source the benchmark file to get the functions
source("benchmark/atime/benchmark-covariance-models-comprehensive.R")

cat("=== TESTING BENCHMARK FRAMEWORK FIXES ===\n")

# Test 1: Data generation produces matrices
cat("Testing data generation...\n")
data_1d <- generate_benchmark_data(20, 1)
data_2d <- generate_benchmark_data(20, 2)

cat(sprintf("1D data dimensions: %s\n", paste(dim(data_1d), collapse = "x")))
cat(sprintf("2D data dimensions: %s\n", paste(dim(data_2d), collapse = "x")))

# Test 2: Multivariate models work
cat("\nTesting multivariate models...\n")
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
test_data <- matrix(rnorm(40), ncol = 2)

for (model in models) {
  tryCatch({
    prior_params <- create_prior_parameters(2, model)
    result <- collect_performance_metrics(model, test_data, prior_params, 10, 5)
    
    if (result$success) {
      cat(sprintf("  ✓ %s: SUCCESS\n", model))
    } else {
      cat(sprintf("  ✗ %s: FAILED - %s\n", model, result$error_message))
    }
  }, error = function(e) {
    cat(sprintf("  ✗ %s: ERROR - %s\n", model, e$message))
  })
}

# Test 3: rWishart wrapper
cat("\nTesting safe rWishart wrapper...\n")
tryCatch({
  # Test with well-conditioned matrix
  Lambda <- diag(2)
  result <- safe_rWishart(1, 3, Lambda)
  cat("  ✓ Well-conditioned matrix: SUCCESS\n")
  
  # Test with ill-conditioned matrix
  Lambda_ill <- matrix(c(1, 0.999, 0.999, 1), 2, 2)
  result <- safe_rWishart(1, 3, Lambda_ill)
  cat("  ✓ Ill-conditioned matrix: SUCCESS\n")
  
}, error = function(e) {
  cat(sprintf("  ✗ Safe rWishart failed: %s\n", e$message))
})

cat("\n=== BENCHMARK FRAMEWORK FIXES TEST COMPLETE ===\n")