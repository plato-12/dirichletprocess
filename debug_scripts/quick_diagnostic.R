# Quick diagnostic to understand the test failure
library(dirichletprocess)
library(testthat)

# Source tolerances
source("tests/testthat/helper-testing.R")

cat("=== Quick Diagnostic ===\n")
cat("CLUSTER_TOLERANCE:", CLUSTER_TOLERANCE, "\n")

# The reported failure shows difference: 0.076 but test fails
# This suggests the actual difference might be 3.076, not 0.076

# Let's simulate what might be happening
reported_diff <- 0.076
actual_diff_possibility <- 3.076

cat("Reported difference:", reported_diff, "\n")
cat("Tolerance:", CLUSTER_TOLERANCE, "\n")
cat("Should pass (reported < tolerance):", reported_diff < CLUSTER_TOLERANCE, "\n")

cat("\nPossible actual difference:", actual_diff_possibility, "\n")
cat("Would fail (actual >= tolerance):", actual_diff_possibility >= CLUSTER_TOLERANCE, "\n")

# Let's check the test failure message format
# "Difference: 0.076" could be the formatted output that truncates digits

# Test a quick consistency run with shorter parameters
cat("\n=== Quick Test Run ===\n")
set.seed(123)
quick_data <- rnorm(20)  # Small sample
quick_results <- validate_r_cpp_consistency("normal", quick_data, iterations = 5, n_runs = 1)

cat("Quick test cluster_count_diff:", quick_results$cluster_count_diff, "\n")
cat("Quick test passes cluster tolerance:", quick_results$cluster_count_diff < CLUSTER_TOLERANCE, "\n")

# Check the rounding/formatting issue
formatted_diff <- sprintf("%.3f", quick_results$cluster_count_diff)
cat("Formatted to 3 decimal places:", formatted_diff, "\n")