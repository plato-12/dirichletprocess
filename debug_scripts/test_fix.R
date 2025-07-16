# Test script to verify the fix for E model initialization

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test function to verify the fix
test_fix <- function() {
  cat("=== Testing Fix for E Model Initialization ===\n")
  
  # Create simple 1D data
  test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  
  # Create E model distribution
  e_model <- MvnormalCreate(list(
    mu0 = c(0),
    kappa0 = 1,
    nu = 2,
    Lambda = diag(1),
    covModel = "E"
  ))
  
  cat("Test data dimensions:", dim(test_data), "\n")
  cat("E model created successfully.\n")
  
  # Test 1: DirichletProcessMvnormal with mixing distribution object (should now work)
  cat("\n1. Testing DirichletProcessMvnormal with mixing distribution object...\n")
  
  tryCatch({
    result1 <- DirichletProcessMvnormal(test_data, e_model)
    cat("SUCCESS: DirichletProcessMvnormal with mixing distribution object works!\n")
    cat("Result clusters:", result1$numberClusters, "\n")
    cat("Result mu dimensions:", dim(result1$clusterParameters$mu), "\n")
    cat("Result sig dimensions:", dim(result1$clusterParameters$sig), "\n")
    
    # Test running a few MCMC steps
    cat("  Testing MCMC steps...\n")
    result1_mcmc <- Fit(result1, its = 10)
    cat("  MCMC successful with", result1_mcmc$iterations, "iterations.\n")
    
  }, error = function(e) {
    cat("FAILED: DirichletProcessMvnormal with mixing distribution object failed:", e$message, "\n")
  })
  
  # Test 2: V model for comparison
  cat("\n2. Testing V model with same approach...\n")
  
  tryCatch({
    v_model <- MvnormalCreate(list(
      mu0 = c(0),
      kappa0 = 1,
      nu = 2,
      Lambda = diag(1),
      covModel = "V"
    ))
    
    result2 <- DirichletProcessMvnormal(test_data, v_model)
    cat("SUCCESS: V model with mixing distribution object works!\n")
    cat("Result clusters:", result2$numberClusters, "\n")
    
  }, error = function(e) {
    cat("FAILED: V model with mixing distribution object failed:", e$message, "\n")
  })
  
  # Test 3: Multivariate models to ensure no regression
  cat("\n3. Testing multivariate models (regression test)...\n")
  
  tryCatch({
    # 2D data
    test_data_2d <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
    
    # EII model
    eii_model <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = "EII"
    ))
    
    result3 <- DirichletProcessMvnormal(test_data_2d, eii_model)
    cat("SUCCESS: EII model works!\n")
    cat("Result clusters:", result3$numberClusters, "\n")
    
  }, error = function(e) {
    cat("FAILED: EII model failed:", e$message, "\n")
  })
  
  cat("\n=== Fix Test Complete ===\n")
}

# Run the test
test_fix()