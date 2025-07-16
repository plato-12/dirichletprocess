# Focus on E model for univariate data (d=1) - the correct usage
# This script analyzes the E model parameter format issues with proper 1D data

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test function to analyze univariate E model
test_univariate_e_model <- function() {
  cat("=== Univariate E Model Analysis ===\n")
  
  # Create simple 1D data (E model requires d=1)
  test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  
  # Create E model distribution for 1D data
  e_model <- MvnormalCreate(list(
    mu0 = c(0),     # 1D mean
    kappa0 = 1,
    nu = 2,         # Minimum nu for d=1
    Lambda = diag(1),  # 1x1 matrix
    covModel = "E"
  ))
  
  # Create V model distribution for 1D data  
  v_model <- MvnormalCreate(list(
    mu0 = c(0),     # 1D mean
    kappa0 = 1,
    nu = 2,         # Minimum nu for d=1
    Lambda = diag(1),  # 1x1 matrix
    covModel = "V"
  ))
  
  cat("Models created successfully for 1D data.\n")
  
  # Test PriorDraw for both models with careful error handling
  cat("\n1. Testing PriorDraw formats...\n")
  
  # Test E model
  cat("Testing E model PriorDraw...\n")
  tryCatch({
    prior_e <- PriorDraw(e_model, n = 3)
    cat("E model PriorDraw successful.\n")
    cat("E model mu class:", class(prior_e$mu), "dim:", dim(prior_e$mu), "\n")
    cat("E model sig class:", class(prior_e$sig), "dim:", dim(prior_e$sig), "\n")
    
    # Check structure
    cat("E model mu structure:\n")
    str(prior_e$mu)
    cat("E model sig structure:\n")
    str(prior_e$sig)
    
  }, error = function(e) {
    cat("E model PriorDraw failed:", e$message, "\n")
  })
  
  # Test V model
  cat("Testing V model PriorDraw...\n")
  tryCatch({
    prior_v <- PriorDraw(v_model, n = 3)
    cat("V model PriorDraw successful.\n")
    cat("V model mu class:", class(prior_v$mu), "dim:", dim(prior_v$mu), "\n")
    cat("V model sig class:", class(prior_v$sig), "dim:", dim(prior_v$sig), "\n")
    
    # Check structure
    cat("V model mu structure:\n")
    str(prior_v$mu)
    cat("V model sig structure:\n")
    str(prior_v$sig)
    
  }, error = function(e) {
    cat("V model PriorDraw failed:", e$message, "\n")
  })
  
  # Test DirichletProcess creation for E model
  cat("\n2. Testing DirichletProcess creation for E model...\n")
  
  tryCatch({
    dp_e <- DirichletProcessMvnormal(test_data, e_model)
    cat("E model DP object created successfully.\n")
    cat("E model initial cluster parameters:\n")
    cat("  mu dimensions:", dim(dp_e$clusterParameters$mu), "\n")
    cat("  sig dimensions:", dim(dp_e$clusterParameters$sig), "\n")
    
    # Test initialization
    cat("Testing E model initialization...\n")
    dp_e_init <- Initialise(dp_e)
    cat("E model initialization successful.\n")
    cat("After initialization - mu dimensions:", dim(dp_e_init$clusterParameters$mu), "\n")
    cat("After initialization - sig dimensions:", dim(dp_e_init$clusterParameters$sig), "\n")
    
  }, error = function(e) {
    cat("E model DP creation/initialization failed:", e$message, "\n")
    cat("Error details:\n")
    print(e)
  })
  
  # Test DirichletProcess creation for V model  
  cat("\n3. Testing DirichletProcess creation for V model...\n")
  
  tryCatch({
    dp_v <- DirichletProcessMvnormal(test_data, v_model)
    cat("V model DP object created successfully.\n")
    cat("V model initial cluster parameters:\n")
    cat("  mu dimensions:", dim(dp_v$clusterParameters$mu), "\n")
    cat("  sig dimensions:", dim(dp_v$clusterParameters$sig), "\n")
    
    # Test initialization
    cat("Testing V model initialization...\n")
    dp_v_init <- Initialise(dp_v)
    cat("V model initialization successful.\n")
    cat("After initialization - mu dimensions:", dim(dp_v_init$clusterParameters$mu), "\n")
    cat("After initialization - sig dimensions:", dim(dp_v_init$clusterParameters$sig), "\n")
    
  }, error = function(e) {
    cat("V model DP creation/initialization failed:", e$message, "\n")
    cat("Error details:\n")
    print(e)
  })
  
  # Test parameter helper functions
  cat("\n4. Testing parameter helper functions...\n")
  
  tryCatch({
    d <- 1
    e_params <- getNumCovParams(d, "E")
    v_params <- getNumCovParams(d, "V")
    
    cat("E model parameters count (d=1):", e_params, "\n")
    cat("V model parameters count (d=1):", v_params, "\n")
    
  }, error = function(e) {
    cat("getNumCovParams failed:", e$message, "\n")
  })
  
  cat("\n=== Analysis Complete ===\n")
}

# Run the test
test_univariate_e_model()