#!/usr/bin/env Rscript

# Test script to verify covariance models implementation
# This tests the basic functionality of all covariance models

library(devtools)
load_all()

# Test function to create and initialize a covariance model
test_covariance_model <- function(model_name, data_dim = 2) {
  cat(sprintf("Testing covariance model: %s\n", model_name))
  
  # Create test data
  if (model_name %in% c("E", "V")) {
    # These are for univariate data only
    data_dim <- 1
  }
  
  set.seed(42)
  if (data_dim == 1) {
    test_data <- rnorm(20)
  } else {
    test_data <- matrix(rnorm(20 * data_dim), ncol = data_dim)
  }
  
  tryCatch({
    # Create mixing distribution
    prior_params <- list(
      mu0 = if (data_dim == 1) 0 else rep(0, data_dim),
      kappa0 = 1,
      nu = data_dim + 1,
      Lambda = if (data_dim == 1) 1 else diag(data_dim),
      covModel = model_name
    )
    md <- MvnormalCreate(prior_params)
    
    # Check class structure
    expected_class <- if (model_name == "FULL") {
      c("list", "mvnormal", "conjugate")
    } else {
      c("list", paste0("mvnormal.", model_name), "mvnormal", "conjugate")
    }
    
    if (!identical(class(md), expected_class)) {
      stop(sprintf("Incorrect class structure. Expected: %s, Got: %s", 
                   paste(expected_class, collapse = ", "), 
                   paste(class(md), collapse = ", ")))
    }
    
    # Create Dirichlet process using DirichletProcessCreate 
    dp <- DirichletProcessCreate(test_data, md)
    
    # Test Initialise method
    dp <- Initialise(dp, numInitialClusters = 2)
    
    # Verify initialization worked
    if (dp$numberClusters != 2) {
      stop(sprintf("Initialization failed: expected 2 clusters, got %d", dp$numberClusters))
    }
    
    # Test basic methods
    PriorDraw(md, 1)
    PosteriorDraw(md, test_data, 1)
    
    cat(sprintf("✓ %s: SUCCESS\n", model_name))
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("✗ %s: FAILED - %s\n", model_name, e$message))
    return(FALSE)
  })
}

# Test all covariance models
models <- c("FULL", "E", "V", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
results <- sapply(models, test_covariance_model)

# Summary
cat("\n=== TEST SUMMARY ===\n")
success_count <- sum(results)
total_count <- length(results)
cat(sprintf("Passed: %d/%d tests\n", success_count, total_count))

if (success_count == total_count) {
  cat("✓ ALL TESTS PASSED! Covariance models are working correctly.\n")
} else {
  cat("✗ Some tests failed. Review the errors above.\n")
  failed_models <- names(results[!results])
  cat(sprintf("Failed models: %s\n", paste(failed_models, collapse = ", ")))
}