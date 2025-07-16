# Diagnostic script for E model initialization issue - R only
# This script tests E model initialization with C++ disabled

library(dirichletprocess)

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test function to analyze E model initialization
test_e_model_r_only <- function() {
  cat("=== E Model R-Only Initialization Diagnostic ===\n")
  
  # Step 1: Create E model distribution
  cat("\n1. Creating E model distribution...\n")
  
  # Create simple 2D data
  test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  
  # Create E model distribution
  e_model <- MvnormalCreate(list(
    mu0 = c(0, 0),
    kappa0 = 1,
    nu = 3,
    Lambda = diag(2),
    covModel = "E"
  ))
  
  cat("E model created successfully.\n")
  cat("E model class:", class(e_model), "\n")
  cat("E model covModel:", e_model$priorParameters$covModel, "\n")
  
  # Step 2: Test PriorDraw with E model
  cat("\n2. Testing PriorDraw with E model...\n")
  
  tryCatch({
    prior_result <- PriorDraw(e_model, n = 3)
    cat("PriorDraw successful.\n")
    cat("mu dimensions:", dim(prior_result$mu), "\n")
    cat("sig dimensions:", dim(prior_result$sig), "\n")
    cat("mu class:", class(prior_result$mu), "\n")
    cat("sig class:", class(prior_result$sig), "\n")
    
    # Show first values
    cat("mu[1,1,1]:", prior_result$mu[1,1,1], "\n")
    cat("sig[1,1]:", prior_result$sig[1,1], "\n")
    
  }, error = function(e) {
    cat("PriorDraw failed with error:", e$message, "\n")
    return(NULL)
  })
  
  # Step 3: Test DirichletProcess creation
  cat("\n3. Testing DirichletProcess creation...\n")
  
  tryCatch({
    dp_obj <- DirichletProcessMvnormal(test_data, e_model)
    cat("DirichletProcess creation successful.\n")
    cat("Number of clusters:", dp_obj$numberClusters, "\n")
    cat("Cluster parameters mu dimensions:", dim(dp_obj$clusterParameters$mu), "\n")
    cat("Cluster parameters sig dimensions:", dim(dp_obj$clusterParameters$sig), "\n")
    
  }, error = function(e) {
    cat("DirichletProcess creation failed with error:", e$message, "\n")
    cat("Error details:\n")
    print(e)
    return(NULL)
  })
  
  # Step 4: Test Initialise function directly
  cat("\n4. Testing Initialise function...\n")
  
  tryCatch({
    # Create minimal DP object
    dp_obj <- DirichletProcessMvnormal(test_data, e_model)
    
    # Test initialization
    dp_initialized <- Initialise(dp_obj)
    cat("Initialise successful.\n")
    cat("After initialization - mu dimensions:", dim(dp_initialized$clusterParameters$mu), "\n")
    cat("After initialization - sig dimensions:", dim(dp_initialized$clusterParameters$sig), "\n")
    
  }, error = function(e) {
    cat("Initialise failed with error:", e$message, "\n")
    cat("Error details:\n")
    print(e)
    return(NULL)
  })
  
  # Step 5: Compare with V model (should work)
  cat("\n5. Comparing with V model (control)...\n")
  
  tryCatch({
    v_model <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = "V"
    ))
    
    dp_v <- DirichletProcessMvnormal(test_data, v_model)
    dp_v_init <- Initialise(dp_v)
    
    cat("V model works correctly.\n")
    cat("V model mu dimensions:", dim(dp_v_init$clusterParameters$mu), "\n")
    cat("V model sig dimensions:", dim(dp_v_init$clusterParameters$sig), "\n")
    
    # Compare parameter formats
    cat("\n=== Parameter Format Comparison ===\n")
    
    # Test PriorDraw for both models
    prior_e <- PriorDraw(e_model, n = 1)
    prior_v <- PriorDraw(v_model, n = 1)
    
    cat("E model PriorDraw mu format:", dim(prior_e$mu), "class:", class(prior_e$mu), "\n")
    cat("V model PriorDraw mu format:", dim(prior_v$mu), "class:", class(prior_v$mu), "\n")
    cat("E model PriorDraw sig format:", dim(prior_e$sig), "class:", class(prior_e$sig), "\n")
    cat("V model PriorDraw sig format:", dim(prior_v$sig), "class:", class(prior_v$sig), "\n")
    
  }, error = function(e) {
    cat("V model comparison failed with error:", e$message, "\n")
  })
  
  cat("\n=== Diagnostic Complete ===\n")
}

# Run the diagnostic
test_e_model_r_only()