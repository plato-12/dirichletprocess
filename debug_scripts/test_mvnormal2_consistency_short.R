# Short MVNormal2 consistency test
library(testthat)
library(dirichletprocess)

# Load tolerance values from helper file
source("tests/testthat/helper-testing.R")

# Test with very short runs - use values from helper-testing.R
# ALPHA_TOLERANCE, CLUSTER_TOLERANCE, LIKELIHOOD_CORR_MIN are now loaded

# Test data
set.seed(123)
mu1 <- c(-2, -2)
mu2 <- c(2, 2)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
test_data <- rbind(
  mvtnorm::rmvnorm(15, mu1, sigma),
  mvtnorm::rmvnorm(15, mu2, sigma)
)

cat("=== Short MVNormal2 Consistency Test ===\n")

# Single run comparison with very short iterations
results <- list()
seed <- 12345

for (run in 1:2) {  # Only 2 runs instead of 5
  current_seed <- seed + run - 1
  
  # R implementation
  set.seed(current_seed)
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessMvnormal2(test_data)
  dp_r <- Fit(dp_r, its = 20, progressBar = FALSE)  # Very short
  
  # C++ implementation  
  set.seed(current_seed)
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessMvnormal2(test_data)
  dp_cpp <- Fit(dp_cpp, its = 20, progressBar = FALSE)  # Very short
  
  # Extract statistics
  r_alpha_mean <- mean(dp_r$alphaChain, na.rm = TRUE)
  cpp_alpha_mean <- mean(dp_cpp$alphaChain, na.rm = TRUE)
  r_clusters <- mean(sapply(dp_r$labelsChain, function(x) length(unique(x))))
  cpp_clusters <- mean(sapply(dp_cpp$labelsChain, function(x) length(unique(x))))
  
  # Store results
  results[[run]] <- list(
    alpha_mean_diff = abs(r_alpha_mean - cpp_alpha_mean),
    cluster_count_diff = abs(r_clusters - cpp_clusters),
    r_alpha = r_alpha_mean,
    cpp_alpha = cpp_alpha_mean,
    r_clusters = r_clusters,
    cpp_clusters = cpp_clusters
  )
  
  cat("Run", run, "- R alpha:", round(r_alpha_mean, 3), 
      "C++ alpha:", round(cpp_alpha_mean, 3),
      "diff:", round(results[[run]]$alpha_mean_diff, 3), "\n")
  cat("Run", run, "- R clusters:", round(r_clusters, 1), 
      "C++ clusters:", round(cpp_clusters, 1),
      "diff:", round(results[[run]]$cluster_count_diff, 1), "\n")
}

# Aggregate results
alpha_mean_diff <- mean(sapply(results, `[[`, "alpha_mean_diff"))
cluster_count_diff <- mean(sapply(results, `[[`, "cluster_count_diff"))

cat("\n=== Final Results ===\n")
cat("Mean alpha difference:", round(alpha_mean_diff, 3), 
    "(tolerance:", ALPHA_TOLERANCE, ") - ", 
    ifelse(alpha_mean_diff < ALPHA_TOLERANCE, "PASS", "FAIL"), "\n")
cat("Mean cluster difference:", round(cluster_count_diff, 1), 
    "(tolerance:", CLUSTER_TOLERANCE, ") - ", 
    ifelse(cluster_count_diff < CLUSTER_TOLERANCE, "PASS", "FAIL"), "\n")

# Test with testthat expectations
cat("\n=== Running testthat checks ===\n")
tryCatch({
  expect_lt(alpha_mean_diff, ALPHA_TOLERANCE)
  cat("✓ Alpha tolerance test PASSED\n")
}, error = function(e) {
  cat("✗ Alpha tolerance test FAILED:", e$message, "\n")
})

tryCatch({
  expect_lt(cluster_count_diff, CLUSTER_TOLERANCE)
  cat("✓ Cluster tolerance test PASSED\n")
}, error = function(e) {
  cat("✗ Cluster tolerance test FAILED:", e$message, "\n")
})

cat("\n=== Test completed ===\n")