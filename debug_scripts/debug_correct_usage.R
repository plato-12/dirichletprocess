# Debug script to test correct usage of DirichletProcessMvnormal

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to test correct usage
debug_correct_usage <- function() {
  cat("=== Debugging Correct Usage of DirichletProcessMvnormal ===\n")
  
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
  
  # Method 1: Pass priorParameters directly
  cat("\n1. Testing DirichletProcessMvnormal with priorParameters...\n")
  
  tryCatch({
    result1 <- DirichletProcessMvnormal(test_data, e_model$priorParameters)
    cat("DirichletProcessMvnormal with priorParameters successful.\n")
    cat("Result clusters:", result1$numberClusters, "\n")
    
  }, error = function(e) {
    cat("DirichletProcessMvnormal with priorParameters failed:", e$message, "\n")
  })
  
  # Method 2: Use DirichletProcessCreate directly
  cat("\n2. Testing DirichletProcessCreate + Initialise...\n")
  
  tryCatch({
    dpobj <- DirichletProcessCreate(test_data, e_model, c(2, 4))
    result2 <- Initialise(dpobj, numInitialClusters = 1)
    cat("DirichletProcessCreate + Initialise successful.\n")
    cat("Result clusters:", result2$numberClusters, "\n")
    
  }, error = function(e) {
    cat("DirichletProcessCreate + Initialise failed:", e$message, "\n")
  })
  
  # Method 3: Test V model for comparison
  cat("\n3. Testing V model with same approach...\n")
  
  tryCatch({
    v_model <- MvnormalCreate(list(
      mu0 = c(0),
      kappa0 = 1,
      nu = 2,
      Lambda = diag(1),
      covModel = "V"
    ))
    
    result3 <- DirichletProcessMvnormal(test_data, v_model$priorParameters)
    cat("V model with priorParameters successful.\n")
    cat("Result clusters:", result3$numberClusters, "\n")
    
  }, error = function(e) {
    cat("V model with priorParameters failed:", e$message, "\n")
  })
  
  # Method 4: Try to reproduce the exact error
  cat("\n4. Reproducing the exact error...\n")
  
  tryCatch({
    # This should fail - passing mixing distribution object instead of priorParameters
    result4 <- DirichletProcessMvnormal(test_data, e_model)
    cat("Unexpected success - this should have failed.\n")
    
  }, error = function(e) {
    cat("Expected error reproduced:", e$message, "\n")
    
    # Now let's see what happens inside MvnormalCreate when passed a mixing distribution
    cat("  Testing MvnormalCreate with mixing distribution object...\n")
    
    tryCatch({
      test_mdobj <- MvnormalCreate(e_model)
      cat("  MvnormalCreate succeeded (unexpected).\n")
      
    }, error = function(e2) {
      cat("  MvnormalCreate failed:", e2$message, "\n")
    })
  })
  
  cat("\n=== Debug Complete ===\n")
}

# Run the debug
debug_correct_usage()