# tests/testthat/test-debug-basic.R
# Basic diagnostic tests

context("Basic Diagnostic Tests")

test_that("System check passes", {
  check_result <- check_cpp_system()
  expect_true(check_result)
})

test_that("Test data generation works", {
  test_data <- create_test_data(n = 30, k = 2)
  expect_true(length(test_data$data) == 30)
  expect_true(test_data$true_k == 2)
})

test_that("Safe comparison works", {
  test_data <- create_test_data(n = 20, k = 2)
  results <- safe_compare_implementations(test_data$data, n_iter = 20)

  expect_true(results$r_success)

  if (results$both_available) {
    inspect_dp_object(results$r, "R")
    inspect_dp_object(results$cpp, "C++")

    # Basic sanity checks
    expect_true(results$r$numberClusters >= 1)
    expect_true((results$cpp$numberClusters %||% results$cpp$n_clusters) >= 1)
  }
})

test_that("Diagnose C++ MCMC chains issue", {
  test_data <- create_test_data(n = 15, k = 2)
  results <- safe_compare_implementations(test_data$data, n_iter = 10)

  skip_if(!results$both_available, "C++ not available")

  cat("\n=== CHAIN DIAGNOSIS ===\n")

  # Check what chains each implementation produces
  r_chains <- c("alphaChain", "likelihoodChain", "weightsChain", "labelsChain")
  r_has_chains <- r_chains[r_chains %in% names(results$r)]
  cpp_has_chains <- r_chains[r_chains %in% names(results$cpp)]

  cat("R has chains:", paste(r_has_chains, collapse = ", "), "\n")
  cat("C++ has chains:", paste(cpp_has_chains, collapse = ", "), "\n")

  if (length(cpp_has_chains) == 0) {
    cat("\nCRITICAL: C++ has NO MCMC chains!\n")
    cat("This means the C++ implementation is not running proper MCMC\n")
    cat("or is not returning the chain data to R correctly.\n")
  }

  expect_true(length(r_has_chains) > 0, "R should have MCMC chains")
})
