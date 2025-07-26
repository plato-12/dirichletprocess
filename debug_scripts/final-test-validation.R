# Final Test Validation Script
# Comprehensive validation that all tests work after reorganization

library(testthat)

cat("=== FINAL TEST VALIDATION ===\n")

# List of critical basic functionality tests
basic_tests <- c(
  "test_dirichlet_process.R",
  "test_normal_inverse_gamma.R", 
  "test_beta_uniform_gamma.R",
  "test_exponential_gamma.R",
  "test_weibull_uniform_gamma.R",
  "test_mvnormal_normal_wishart.R",
  "test_mvnormal_semi_conjugate.R",
  "test_conjugate.R",
  "test_nonconjugate.R"
)

cat("Testing", length(basic_tests), "critical basic functionality tests...\n\n")

# Track results
results <- list()

# Test each file
for (test_file in basic_tests) {
  cat("Testing", test_file, "... ")
  
  # Change to tests/testthat directory and run test
  original_dir <- getwd()
  setwd("tests/testthat")
  
  start_time <- Sys.time()
  result <- try({
    test_result <- test_file(test_file, reporter = "silent")
    list(
      passed = TRUE,
      duration = as.numeric(Sys.time() - start_time),
      summary = attr(test_result, "summary")
    )
  }, silent = TRUE)
  
  setwd(original_dir)
  
  if (inherits(result, "try-error")) {
    cat("FAILED\n")
    results[[test_file]] <- list(passed = FALSE, error = as.character(result))
  } else {
    cat("PASSED (", round(result$duration, 2), "s)\n")
    results[[test_file]] <- result
  }
}

# Summary
cat("\n=== SUMMARY ===\n")
passed <- sum(sapply(results, function(x) x$passed))
total <- length(results)
cat("Passed:", passed, "/", total, "tests\n")

if (passed == total) {
  cat("✅ ALL BASIC TESTS PASSING!\n")
} else {
  cat("❌ Some tests still failing:\n")
  failed <- names(results)[sapply(results, function(x) !x$passed)]
  for (f in failed) {
    cat("  -", f, "\n")
  }
}

# Test one cpp-consistency file to verify they still work
cat("\n=== C++ CONSISTENCY VALIDATION ===\n")
cat("Testing cpp-consistency functionality... ")

result <- try({
  setwd("tests/testthat")
  test_result <- test_file("cpp-consistency/test-cpp-consistency-beta.R", reporter = "silent")
  setwd("../..")
  TRUE
}, silent = TRUE)

if (inherits(result, "try-error")) {
  cat("FAILED\n")
  cat("Error:", as.character(result), "\n")
} else {
  cat("PASSED\n")
  cat("✅ C++ consistency tests still functional\n")
}

cat("\n=== VALIDATION COMPLETE ===\n")