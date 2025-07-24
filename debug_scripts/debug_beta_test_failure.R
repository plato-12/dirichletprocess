# Debug script for beta distribution test failure
# 2025-07-24

library(dirichletprocess)
source('tests/testthat/helper-testing.R')

# Turn off progress bars to reduce output
options(dirichletprocess.show_progress = FALSE)

# Function to run a single beta consistency test and capture results
debug_beta_consistency <- function(seed = 123, iterations = 50) {
  set.seed(seed)
  
  # Create test data similar to the failing test
  beta_size <- 25
  test_data <- c(rbeta(beta_size, 2, 5), rbeta(beta_size, 5, 2))
  
  cat("Running beta consistency test with seed:", seed, "\n")
  cat("Data summary:\n")
  print(summary(test_data))
  
  # Run the consistency validation
  results <- validate_r_cpp_consistency("beta", test_data, iterations = iterations)
  
  # Print detailed results
  cat("\n=== RESULTS ===\n")
  cat("Alpha mean diff:", results$alpha_mean_diff, "(tolerance:", ALPHA_TOLERANCE, ")\n")
  cat("Cluster count diff:", results$cluster_count_diff, "(tolerance:", CLUSTER_TOLERANCE, ")\n")
  cat("Likelihood correlation:", results$likelihood_correlation, "(min:", LIKELIHOOD_CORR_MIN, ")\n")
  cat("Parameter max diff:", results$param_max_diff, "(tolerance:", PARAM_TOLERANCE, ")\n")
  
  cat("\n=== TEST RESULTS ===\n")
  cat("Alpha test:", ifelse(results$alpha_mean_diff < ALPHA_TOLERANCE, "PASS", "FAIL"), "\n")
  cat("Cluster test:", ifelse(results$cluster_count_diff < CLUSTER_TOLERANCE, "PASS", "FAIL"), "\n")
  cat("Likelihood test:", ifelse(results$likelihood_correlation > LIKELIHOOD_CORR_MIN, "PASS", "FAIL"), "\n")
  
  # Show individual run results
  cat("\n=== INDIVIDUAL RUN RESULTS ===\n")
  for (i in 1:length(results$all_runs)) {
    run <- results$all_runs[[i]]
    cat("Run", i, "- Cluster diff:", run$cluster_count_diff, "\n")
  }
  
  return(results)
}

# Run the debug
cat("Debugging beta consistency test failure...\n")
results <- debug_beta_consistency()

# Additional investigation: run multiple times to see variability
cat("\n=== RUNNING MULTIPLE TESTS TO CHECK VARIABILITY ===\n")
multiple_results <- list()
for (i in 1:3) {
  cat("\n--- Test", i, "---\n")
  multiple_results[[i]] <- debug_beta_consistency(seed = 123 + i, iterations = 30)
}

# Summary of multiple runs
cluster_diffs <- sapply(multiple_results, function(x) x$cluster_count_diff)
cat("\nCluster differences across multiple runs:", paste(round(cluster_diffs, 3), collapse = ", "), "\n")
cat("Max cluster difference:", max(cluster_diffs), "\n")
cat("Current tolerance:", CLUSTER_TOLERANCE, "\n")
cat("Suggested tolerance:", max(cluster_diffs) * 1.1, "\n")