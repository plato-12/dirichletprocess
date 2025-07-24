# Analyze cluster count variance to determine if tolerance needs adjustment
# D:\dirichletprocess\dirichletprocess\debug_scripts\analyze_cluster_variance.R

library(dirichletprocess)
library(testthat)

# Source helper functions
source("tests/testthat/helper-testing.R")

cat("=== Analyzing Cluster Count Variance Across Multiple Test Runs ===\n")
cat("CURRENT CLUSTER_TOLERANCE:", CLUSTER_TOLERANCE, "\n\n")

# Parameters matching the failing test
BASE_ITERATIONS <- 50
BASE_SAMPLE_SIZE <- 50

# Run multiple test instances with different seeds
test_results <- list()
seeds <- c(123, 456, 789, 1011, 1213, 1415, 1617, 1819, 2021, 2223)

for (i in seq_along(seeds)) {
  cat("Running test", i, "with seed", seeds[i], "...\n")
  
  # Generate test data with this seed
  set.seed(seeds[i])
  test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 0, 2), sd = 1)
  
  # Run consistency test
  tryCatch({
    results <- validate_r_cpp_consistency("normal", test_data, iterations = BASE_ITERATIONS)
    test_results[[i]] <- list(
      seed = seeds[i],
      cluster_count_diff = results$cluster_count_diff,
      alpha_mean_diff = results$alpha_mean_diff,
      success = TRUE
    )
    cat("  Cluster count diff:", results$cluster_count_diff, "\n")
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    test_results[[i]] <- list(
      seed = seeds[i],
      cluster_count_diff = NA,
      alpha_mean_diff = NA,
      success = FALSE,
      error = e$message
    )
  })
}

# Analyze results
successful_tests <- test_results[sapply(test_results, function(x) x$success)]
cat("\n=== Analysis Results ===\n")
cat("Successful tests:", length(successful_tests), "out of", length(test_results), "\n")

if (length(successful_tests) > 0) {
  cluster_diffs <- sapply(successful_tests, function(x) x$cluster_count_diff)
  alpha_diffs <- sapply(successful_tests, function(x) x$alpha_mean_diff)
  
  cat("\nCluster Count Differences:\n")
  cat("  Values:", paste(round(cluster_diffs, 3), collapse = ", "), "\n")
  cat("  Mean:", round(mean(cluster_diffs), 3), "\n")
  cat("  Median:", round(median(cluster_diffs), 3), "\n")
  cat("  SD:", round(sd(cluster_diffs), 3), "\n")
  cat("  Min:", round(min(cluster_diffs), 3), "\n")
  cat("  Max:", round(max(cluster_diffs), 3), "\n")
  cat("  95th percentile:", round(quantile(cluster_diffs, 0.95), 3), "\n")
  cat("  99th percentile:", round(quantile(cluster_diffs, 0.99), 3), "\n")
  
  # Check how many exceed current tolerance
  exceeding_tolerance <- sum(cluster_diffs >= CLUSTER_TOLERANCE)
  cat("\nTests exceeding current tolerance (", CLUSTER_TOLERANCE, "):", exceeding_tolerance, "/", length(cluster_diffs), "\n")
  cat("Failure rate:", round(100 * exceeding_tolerance / length(cluster_diffs), 1), "%\n")
  
  # Suggest new tolerance based on data
  suggested_tolerance_95 <- quantile(cluster_diffs, 0.95)
  suggested_tolerance_99 <- quantile(cluster_diffs, 0.99)
  suggested_tolerance_max <- max(cluster_diffs) * 1.1  # 10% buffer
  
  cat("\nSuggested tolerance values:\n")
  cat("  For 95% pass rate:", round(suggested_tolerance_95, 1), "\n")
  cat("  For 99% pass rate:", round(suggested_tolerance_99, 1), "\n")
  cat("  Based on max + 10% buffer:", round(suggested_tolerance_max, 1), "\n")
  
  # Test the original failing case
  cat("\n=== Testing Original Failing Case (seed 123) ===\n")
  original_case <- successful_tests[[which(sapply(successful_tests, function(x) x$seed) == 123)]]
  if (!is.null(original_case)) {
    cat("Cluster count diff:", original_case$cluster_count_diff, "\n")
    cat("Exceeds tolerance:", original_case$cluster_count_diff >= CLUSTER_TOLERANCE, "\n")
    cat("Would pass with tolerance", round(suggested_tolerance_95, 1), ":", original_case$cluster_count_diff < suggested_tolerance_95, "\n")
  }
} else {
  cat("No successful tests - there may be a systematic issue\n")
}

cat("\n=== RECOMMENDATION ===\n")
if (length(successful_tests) > 0) {
  cluster_diffs <- sapply(successful_tests, function(x) x$cluster_count_diff)
  max_diff <- max(cluster_diffs)
  
  if (max_diff > CLUSTER_TOLERANCE) {
    recommended_tolerance <- ceiling(max_diff * 1.2)  # 20% buffer
    cat("Current tolerance (", CLUSTER_TOLERANCE, ") is too strict.\n")
    cat("Recommend increasing CLUSTER_TOLERANCE to:", recommended_tolerance, "\n")
    cat("This would accommodate the observed variance in MCMC clustering.\n")
  } else {
    cat("Current tolerance appears adequate based on this sample.\n")
    cat("The failure may be due to random variance or a specific seed issue.\n")
  }
} else {
  cat("Unable to provide recommendation due to systematic test failures.\n")
}