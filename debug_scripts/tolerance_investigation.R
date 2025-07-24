# Investigate the tolerance issue
library(dirichletprocess)
library(testthat)

source("tests/testthat/helper-testing.R")

cat("=== Tolerance Investigation ===\n")
cat("CLUSTER_TOLERANCE =", CLUSTER_TOLERANCE, "\n")

# Test the specific failing values
failing_values <- c(1.44, 0.967, 0.45)

for(val in failing_values) {
  cat("\nTesting value:", val, "\n")
  cat("  val < CLUSTER_TOLERANCE:", val < CLUSTER_TOLERANCE, "\n")
  cat("  expect_lt would pass:", val < CLUSTER_TOLERANCE, "\n")
  
  # Test the exact testthat condition
  tryCatch({
    expect_lt(val, CLUSTER_TOLERANCE)
    cat("  expect_lt: PASS\n")
  }, error = function(e) {
    cat("  expect_lt: FAIL -", e$message, "\n")
  })
}

# Check if there's a different CLUSTER_TOLERANCE being used
cat("\n=== Environment Check ===\n")
cat("Current CLUSTER_TOLERANCE:", CLUSTER_TOLERANCE, "\n")
cat("Type:", typeof(CLUSTER_TOLERANCE), "\n")
cat("Class:", class(CLUSTER_TOLERANCE), "\n")

# Check if variable exists in global environment
if(exists("CLUSTER_TOLERANCE", envir = .GlobalEnv)) {
  cat("Global CLUSTER_TOLERANCE:", get("CLUSTER_TOLERANCE", envir = .GlobalEnv), "\n")
} else {
  cat("No global CLUSTER_TOLERANCE found\n")
}

# Let's manually run a small version of the validation
cat("\n=== Manual Validation Test ===\n")
set.seed(123)
test_data <- rnorm(25)
quick_result <- validate_r_cpp_consistency("normal", test_data, iterations = 10, n_runs = 1)
cat("Quick cluster_count_diff:", quick_result$cluster_count_diff, "\n")
cat("Quick test result < tolerance:", quick_result$cluster_count_diff < CLUSTER_TOLERANCE, "\n")