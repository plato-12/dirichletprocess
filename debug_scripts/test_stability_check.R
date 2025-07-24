# Test Stability Check - Run the exact failing test multiple times
library(dirichletprocess)
library(testthat)

source("tests/testthat/helper-testing.R")

cat("=== Test Stability Check ===\n")
cat("Running the exact failing test scenario multiple times...\n")

# Exact test parameters
DEV_MODE <- Sys.getenv("DP_DEV_TESTING", "TRUE") == "TRUE"
BASE_ITERATIONS <- if (DEV_MODE) 50 else 200
BASE_SAMPLE_SIZE <- if (DEV_MODE) 50 else 100

cat("Parameters: DEV_MODE =", DEV_MODE, ", iterations =", BASE_ITERATIONS, ", sample_size =", BASE_SAMPLE_SIZE, "\n")

# Run the test 3 times with the same seed to see consistency
results_list <- list()

for(i in 1:3) {
  cat("\n--- Run", i, "---\n")
  
  # Use the exact same seed as the test
  set.seed(123)
  test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 0, 2), sd = 1)
  
  # Run the consistency test
  results <- validate_r_cpp_consistency("normal", test_data, iterations = BASE_ITERATIONS)
  results_list[[i]] <- results
  
  cat("cluster_count_diff:", results$cluster_count_diff, "\n")
  cat("Passes test:", results$cluster_count_diff < CLUSTER_TOLERANCE, "\n")
  cat("alpha_mean_diff:", results$alpha_mean_diff, "\n")
  cat("likelihood_correlation:", results$likelihood_correlation, "\n")
}

# Check if results are consistent across runs
cat("\n=== Consistency Analysis ===\n")
cluster_diffs <- sapply(results_list, function(x) x$cluster_count_diff)
cat("Cluster differences across runs:", cluster_diffs, "\n")
cat("All the same?", length(unique(round(cluster_diffs, 6))) == 1, "\n")

# Check if any run fails
failing_runs <- which(cluster_diffs >= CLUSTER_TOLERANCE)
cat("Failing runs:", if(length(failing_runs) == 0) "None" else failing_runs, "\n")
cat("Failure rate:", length(failing_runs), "/", length(results_list), "\n")