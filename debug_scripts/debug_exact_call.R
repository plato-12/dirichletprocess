# Debug script to make the exact same call as DirichletProcessMvnormal

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to make exact same call
debug_exact_call <- function() {
  cat("=== Debugging Exact DirichletProcessMvnormal Call ===\n")
  
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
  
  # Call DirichletProcessMvnormal exactly as it's called in our tests
  cat("\n1. Testing DirichletProcessMvnormal with exact call...\n")
  
  tryCatch({
    result <- DirichletProcessMvnormal(test_data, e_model)
    cat("DirichletProcessMvnormal successful.\n")
    
  }, error = function(e) {
    cat("DirichletProcessMvnormal failed:", e$message, "\n")
    print(e)
  })
  
  # Now try passing g0Priors directly as the function expects
  cat("\n2. Testing DirichletProcessMvnormal with g0Priors...\n")
  
  tryCatch({
    result2 <- DirichletProcessMvnormal(test_data, list(
      mu0 = c(0),
      kappa0 = 1,
      nu = 2,
      Lambda = diag(1),
      covModel = "E"
    ))
    cat("DirichletProcessMvnormal with g0Priors successful.\n")
    
  }, error = function(e) {
    cat("DirichletProcessMvnormal with g0Priors failed:", e$message, "\n")
    print(e)
  })
  
  # Check the actual source code of DirichletProcessMvnormal
  cat("\n3. Checking DirichletProcessMvnormal source...\n")
  
  tryCatch({
    # Get the source to see what's happening
    cat("DirichletProcessMvnormal function:\n")
    print(DirichletProcessMvnormal)
    
  }, error = function(e) {
    cat("Could not get DirichletProcessMvnormal source:", e$message, "\n")
  })
  
  cat("\n=== Debug Complete ===\n")
}

# Run the debug
debug_exact_call()