# Quick test to verify beta consistency with new tolerance
library(dirichletprocess)
source('tests/testthat/helper-testing.R')

# Disable progress bars
options(dirichletprocess.show_progress = FALSE)

set.seed(123)
beta_size <- 25
test_data <- c(rbeta(beta_size, 2, 5), rbeta(beta_size, 5, 2))

cat("Running beta consistency test...\n")
results <- validate_r_cpp_consistency("beta", test_data, iterations = 50)

cat("Results:\n")
cat("- Alpha mean diff:", results$alpha_mean_diff, "(tolerance:", ALPHA_TOLERANCE, ")\n")
cat("- Cluster count diff:", results$cluster_count_diff, "(tolerance:", CLUSTER_TOLERANCE, ")\n") 
cat("- Likelihood correlation:", results$likelihood_correlation, "(min:", LIKELIHOOD_CORR_MIN, ")\n")

cat("\nTest outcomes:\n")
cat("- Alpha test:", ifelse(results$alpha_mean_diff < ALPHA_TOLERANCE, "PASS", "FAIL"), "\n")
cat("- Cluster test:", ifelse(results$cluster_count_diff < CLUSTER_TOLERANCE, "PASS", "FAIL"), "\n")
cat("- Likelihood test:", ifelse(results$likelihood_correlation > LIKELIHOOD_CORR_MIN, "PASS", "FAIL"), "\n")

cat("\nOverall:", ifelse(
  results$alpha_mean_diff < ALPHA_TOLERANCE &&
  results$cluster_count_diff < CLUSTER_TOLERANCE &&
  results$likelihood_correlation > LIKELIHOOD_CORR_MIN,
  "ALL TESTS PASS", "SOME TESTS FAIL"
), "\n")