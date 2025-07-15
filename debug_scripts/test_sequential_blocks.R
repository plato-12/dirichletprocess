# Test Sequential Execution of Multiple Blocks
# This tests if crashes occur when running multiple test blocks together

library(dirichletprocess)
library(testthat)
library(mvtnorm)

set_use_cpp(TRUE)

cat("=== Testing Sequential Block Execution ===\n")

# Source the comprehensive test functions
source("tests/testthat/test-mvnormal-cpp-comprehensive.R", local = TRUE)

cat("Comprehensive test file sourced successfully\n")

# Try running testthat on the file directly
cat("Attempting to run testthat on the comprehensive file...\n")

tryCatch({
  # This should replicate what testthat::test_file() does
  result <- testthat::test_file("tests/testthat/test-mvnormal-cpp-comprehensive.R")
  cat("✓ testthat::test_file() completed successfully\n")
  print(result)
}, error = function(e) {
  cat("✗ testthat::test_file() failed:", e$message, "\n")
})

set_use_cpp(FALSE)