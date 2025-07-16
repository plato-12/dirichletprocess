# Debug script to trace DirichletProcessCreate issue with E/V models

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to trace DirichletProcessCreate
debug_dirichletprocess_create <- function() {
  cat("=== Debugging DirichletProcessCreate for E/V Models ===\n")
  
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
  
  cat("E model created successfully.\n")
  
  # Step 1: Test DirichletProcessCreate directly
  cat("\n1. Testing DirichletProcessCreate directly...\n")
  
  tryCatch({
    dpobj <- DirichletProcessCreate(test_data, e_model, c(2, 4))
    cat("DirichletProcessCreate successful.\n")
    cat("DP object numberClusters:", dpobj$numberClusters, "\n")
    
    # Check cluster parameters
    cat("Cluster parameters mu dimensions:", dim(dpobj$clusterParameters$mu), "\n")
    cat("Cluster parameters sig dimensions:", dim(dpobj$clusterParameters$sig), "\n")
    
  }, error = function(e) {
    cat("DirichletProcessCreate failed:", e$message, "\n")
    cat("Error details:\n")
    print(e)
    
    # Let's manually trace DirichletProcessCreate
    cat("\nManual trace of DirichletProcessCreate:\n")
    
    # Create base object
    tryCatch({
      base_dp <- list(
        data = test_data,
        mixingDistribution = e_model,
        alphaPriorParameters = c(2, 4),
        numInitialClusters = 1
      )
      
      cat("Base DP object created.\n")
      
      # Test creating empty cluster parameters
      empty_mu <- array(NA, dim = c(1, 1, 1))
      empty_sig <- array(NA, dim = c(1, 1, 1))
      
      cat("Empty cluster parameters created.\n")
      
      # Test PosteriorDraw
      post_result <- PosteriorDraw(e_model, test_data[1:2, , drop = FALSE], n = 1)
      cat("PosteriorDraw successful.\n")
      cat("PosteriorDraw mu dimensions:", dim(post_result$mu), "\n")
      cat("PosteriorDraw sig dimensions:", dim(post_result$sig), "\n")
      
    }, error = function(e2) {
      cat("Manual trace failed:", e2$message, "\n")
    })
  })
  
  # Step 2: Test with DirichletProcessMvnormal (the wrapper)
  cat("\n2. Testing DirichletProcessMvnormal wrapper...\n")
  
  tryCatch({
    dpobj2 <- DirichletProcessMvnormal(test_data, e_model)
    cat("DirichletProcessMvnormal successful.\n")
    
  }, error = function(e) {
    cat("DirichletProcessMvnormal failed:", e$message, "\n")
    cat("Error details:\n")
    print(e)
    traceback()
  })
  
  cat("\n=== Debug Complete ===\n")
}

# Run the debug
debug_dirichletprocess_create()