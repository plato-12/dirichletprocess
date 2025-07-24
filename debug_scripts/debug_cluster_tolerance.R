# Debug script for cluster tolerance issue
# D:\dirichletprocess\dirichletprocess\debug_scripts\debug_cluster_tolerance.R

library(dirichletprocess)
library(testthat)

# Source helper functions
source("tests/testthat/helper-testing.R")

# Set development mode
Sys.setenv(DP_DEV_TESTING = "TRUE")
DEV_MODE <- TRUE
BASE_ITERATIONS <- 50
BASE_SAMPLE_SIZE <- 50

cat("=== Debugging Cluster Tolerance Issue ===\n")
cat("CLUSTER_TOLERANCE:", CLUSTER_TOLERANCE, "\n")
cat("BASE_ITERATIONS:", BASE_ITERATIONS, "\n")
cat("BASE_SAMPLE_SIZE:", BASE_SAMPLE_SIZE, "\n\n")

# Generate test data (same as in the test)
set.seed(123)
test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 0, 2), sd = 1)

cat("Test data summary:\n")
print(summary(test_data))
cat("\n")

# Run consistency tests with detailed output
cat("Running consistency test...\n")
results <- validate_r_cpp_consistency("normal", test_data, iterations = BASE_ITERATIONS)

cat("\n=== Results ===\n")
cat("Alpha mean diff:", results$alpha_mean_diff, "(tolerance:", ALPHA_TOLERANCE, ")\n")
cat("Cluster count diff:", results$cluster_count_diff, "(tolerance:", CLUSTER_TOLERANCE, ")\n")
cat("Likelihood correlation:", results$likelihood_correlation, "(min:", LIKELIHOOD_CORR_MIN, ")\n")
cat("Param max diff:", results$param_max_diff, "(tolerance:", PARAM_TOLERANCE, ")\n")

cat("\n=== Test Status ===\n")
cat("Alpha test:", if(results$alpha_mean_diff < ALPHA_TOLERANCE) "PASS" else "FAIL", "\n")
cat("Cluster test:", if(results$cluster_count_diff < CLUSTER_TOLERANCE) "PASS" else "FAIL", "\n")
cat("Likelihood test:", if(results$likelihood_correlation > LIKELIHOOD_CORR_MIN) "PASS" else "FAIL", "\n")
cat("Param test:", if(results$param_max_diff < PARAM_TOLERANCE) "PASS" else "FAIL", "\n")

# Look at individual run results
cat("\n=== Individual Run Results ===\n")
for (i in seq_along(results$all_runs)) {
  run <- results$all_runs[[i]]
  cat("Run", i, ": cluster_count_diff =", run$cluster_count_diff, "\n")
}

cat("\n=== Analysis ===\n")
cluster_diffs <- sapply(results$all_runs, function(x) x$cluster_count_diff)
cat("Cluster count differences across runs:", paste(cluster_diffs, collapse = ", "), "\n")
cat("Mean cluster diff:", mean(cluster_diffs), "\n")
cat("Max cluster diff:", max(cluster_diffs), "\n")
cat("Min cluster diff:", min(cluster_diffs), "\n")

# Check if it's a systematic issue or just variance
if (results$cluster_count_diff >= CLUSTER_TOLERANCE) {
  cat("\nFAILURE ANALYSIS:\n")
  cat("The cluster count difference (", results$cluster_count_diff, ") exceeds tolerance (", CLUSTER_TOLERANCE, ")\n")
  cat("Difference margin:", results$cluster_count_diff - CLUSTER_TOLERANCE, "\n")
  
  if (max(cluster_diffs) > CLUSTER_TOLERANCE) {
    cat("Issue appears to be in individual runs, not just aggregation\n")
  } else {
    cat("Issue appears to be in aggregation across runs\n")
  }
}