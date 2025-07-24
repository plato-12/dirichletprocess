# Debug the specific test failure
library(dirichletprocess)
library(testthat)

# Source the helper functions to get tolerances
source("tests/testthat/helper-testing.R")

cat("=== Debug Test Failure ===\n")
cat("CLUSTER_TOLERANCE:", CLUSTER_TOLERANCE, "\n")

# Recreate the exact test scenario
DEV_MODE <- Sys.getenv("DP_DEV_TESTING", "TRUE") == "TRUE"
BASE_ITERATIONS <- if (DEV_MODE) 50 else 200
BASE_SAMPLE_SIZE <- if (DEV_MODE) 50 else 100

cat("DEV_MODE:", DEV_MODE, "\n")
cat("BASE_ITERATIONS:", BASE_ITERATIONS, "\n")
cat("BASE_SAMPLE_SIZE:", BASE_SAMPLE_SIZE, "\n")

# Generate test data exactly as in the test
set.seed(123)
test_data <- rnorm(BASE_SAMPLE_SIZE, mean = c(-2, 0, 2), sd = 1)

# Run consistency test
results <- validate_r_cpp_consistency("normal", test_data, iterations = BASE_ITERATIONS)

cat("\n=== Results ===\n")
cat("alpha_mean_diff:", results$alpha_mean_diff, "(tolerance:", ALPHA_TOLERANCE, ")\n")
cat("cluster_count_diff:", results$cluster_count_diff, "(tolerance:", CLUSTER_TOLERANCE, ")\n")
cat("likelihood_correlation:", results$likelihood_correlation, "(min:", LIKELIHOOD_CORR_MIN, ")\n")
cat("param_max_diff:", results$param_max_diff, "(tolerance:", PARAM_TOLERANCE, ")\n")

# Check each condition
cat("\n=== Test Conditions ===\n")
cat("alpha_mean_diff < ALPHA_TOLERANCE:", results$alpha_mean_diff < ALPHA_TOLERANCE, "\n")
cat("cluster_count_diff < CLUSTER_TOLERANCE:", results$cluster_count_diff < CLUSTER_TOLERANCE, "\n")
cat("likelihood_correlation > LIKELIHOOD_CORR_MIN:", results$likelihood_correlation > LIKELIHOOD_CORR_MIN, "\n")
cat("param_max_diff < PARAM_TOLERANCE:", results$param_max_diff < PARAM_TOLERANCE, "\n")

# Check if the issue is with expect_lt vs expect_lte
cat("\n=== Detailed Analysis ===\n")
cat("Cluster difference:", results$cluster_count_diff, "\n")
cat("Is exactly equal to tolerance?", results$cluster_count_diff == CLUSTER_TOLERANCE, "\n")
cat("Difference from tolerance:", results$cluster_count_diff - CLUSTER_TOLERANCE, "\n")

# Test the exact condition used in testthat
cat("Using expect_lt logic (strictly less than):", results$cluster_count_diff < CLUSTER_TOLERANCE, "\n")