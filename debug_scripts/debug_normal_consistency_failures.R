# Debug script for Normal distribution consistency test failures
# This script examines the actual values that are causing test failures

library(dirichletprocess)
library(testthat)

# Source the helper functions
source("tests/testthat/helper-testing.R")

# Set development mode for faster debugging
DEV_MODE <- TRUE
BASE_ITERATIONS <- 50
BASE_SAMPLE_SIZE <- 50

cat("=== DEBUG: Normal Distribution Consistency Test Failures ===\n")
cat("CLUSTER_TOLERANCE:", CLUSTER_TOLERANCE, "\n")
cat("ALPHA_TOLERANCE:", ALPHA_TOLERANCE, "\n")
cat("Base iterations:", BASE_ITERATIONS, "\n")
cat("Base sample size:", BASE_SAMPLE_SIZE, "\n\n")

# Test 1: Basic consistency test
cat("--- Test 1: Basic Consistency ---\n")
set.seed(123)
test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 0, 2), sd = 1)
cat("Test data length:", length(test_data), "\n")
cat("Test data summary:", summary(test_data), "\n")

results <- validate_r_cpp_consistency("normal", test_data, iterations = BASE_ITERATIONS)

cat("Results:\n")
cat("  alpha_mean_diff:", results$alpha_mean_diff, "(tolerance:", ALPHA_TOLERANCE, ")\n")
cat("  cluster_count_diff:", results$cluster_count_diff, "(tolerance:", CLUSTER_TOLERANCE, ")\n")
cat("  likelihood_correlation:", results$likelihood_correlation, "(min:", LIKELIHOOD_CORR_MIN, ")\n")
cat("  param_max_diff:", results$param_max_diff, "(tolerance:", PARAM_TOLERANCE, ")\n")

cat("\nTest 1 Failures:\n")
if (results$cluster_count_diff >= CLUSTER_TOLERANCE) {
  cat("  FAIL: cluster_count_diff =", results$cluster_count_diff, ">= tolerance", CLUSTER_TOLERANCE, "\n")
}
if (results$alpha_mean_diff >= ALPHA_TOLERANCE) {
  cat("  FAIL: alpha_mean_diff =", results$alpha_mean_diff, ">= tolerance", ALPHA_TOLERANCE, "\n")
}

# Test 2: Different sample sizes
cat("\n--- Test 2: Different Sample Sizes ---\n")
sample_sizes <- c(25, 50)
dev_iterations <- 30
dev_runs <- 2

for (n in sample_sizes) {
  cat("\nTesting sample size:", n, "\n")
  test_data <- generate_test_data("normal", n)
  cat("Test data length:", length(test_data), "\n")
  
  results <- validate_r_cpp_consistency("normal", test_data, 
                                       iterations = dev_iterations, 
                                       n_runs = dev_runs)
  
  cat("Results for n =", n, ":\n")
  cat("  alpha_mean_diff:", results$alpha_mean_diff, "(tolerance:", ALPHA_TOLERANCE, ")\n")
  cat("  cluster_count_diff:", results$cluster_count_diff, "(tolerance:", CLUSTER_TOLERANCE, ")\n")
  
  if (results$cluster_count_diff >= CLUSTER_TOLERANCE) {
    cat("  FAIL: cluster_count_diff =", results$cluster_count_diff, ">= tolerance", CLUSTER_TOLERANCE, "\n")
  }
  if (results$alpha_mean_diff >= ALPHA_TOLERANCE) {
    cat("  FAIL: alpha_mean_diff =", results$alpha_mean_diff, ">= tolerance", ALPHA_TOLERANCE, "\n")
  }
}

# Test 3: Individual run analysis for n=25
cat("\n--- Test 3: Individual Run Analysis (n=25) ---\n")
test_data <- generate_test_data("normal", 25)
results <- validate_r_cpp_consistency("normal", test_data, iterations = 30, n_runs = 2)

cat("Individual run details:\n")
for (i in 1:length(results$all_runs)) {
  run_result <- results$all_runs[[i]]
  cat("Run", i, ":\n")
  cat("  alpha_mean_diff:", run_result$alpha_mean_diff, "\n")
  cat("  cluster_count_diff:", run_result$cluster_count_diff, "\n")
  cat("  likelihood_correlation:", run_result$likelihood_correlation, "\n")
}

cat("\n=== SUMMARY ===\n")
cat("The failures appear to be due to:\n")
cat("1. cluster_count_diff exceeding tolerance of", CLUSTER_TOLERANCE, "\n")
cat("2. This suggests R and C++ implementations are producing different cluster structures\n")
cat("3. The issue may be:\n")
cat("   - Different random number generation between R and C++\n")
cat("   - Different algorithmic choices in cluster assignment\n")
cat("   - Initialization differences\n")
cat("   - Parameter handling differences\n")