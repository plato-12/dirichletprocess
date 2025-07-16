# Focus on parameter format differences between E and V models - v2
# This script analyzes the exact parameter format issues with proper loading

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test function to analyze parameter formats
test_parameter_formats <- function() {
  cat("=== Parameter Format Analysis v2 ===\n")
  
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
  
  # Create V model distribution  
  v_model <- MvnormalCreate(list(
    mu0 = c(0, 0),
    kappa0 = 1,
    nu = 3,
    Lambda = diag(2),
    covModel = "V"
  ))
  
  cat("Models created successfully.\n")
  
  # Test PriorDraw for both models with careful error handling
  cat("\n1. Testing PriorDraw formats...\n")
  
  # Test E model
  cat("Testing E model PriorDraw...\n")
  tryCatch({
    prior_e <- PriorDraw(e_model, n = 1)
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
    prior_v <- PriorDraw(v_model, n = 1)
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
  
  # Test parameter helper functions
  cat("\n2. Testing parameter helper functions...\n")
  
  tryCatch({
    d <- 2
    e_params <- getNumCovParams(d, "E")
    v_params <- getNumCovParams(d, "V")
    
    cat("E model parameters count:", e_params, "\n")
    cat("V model parameters count:", v_params, "\n")
    
  }, error = function(e) {
    cat("getNumCovParams failed:", e$message, "\n")
  })
  
  # Test parameter extraction/reconstruction  
  cat("\n3. Testing parameter extraction/reconstruction...\n")
  
  tryCatch({
    # Create a simple covariance matrix
    test_cov <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
    
    e_extracted <- extractCovarianceParams(test_cov, "E")
    v_extracted <- extractCovarianceParams(test_cov, "V")
    
    cat("E model extracted params:", e_extracted, "\n")
    cat("V model extracted params:", v_extracted, "\n")
    
    # Test reconstruction
    e_reconstructed <- reconstructCovarianceMatrix(e_extracted, 2, "E")
    v_reconstructed <- reconstructCovarianceMatrix(v_extracted, 2, "V")
    
    cat("E model reconstruction successful.\n")
    cat("V model reconstruction successful.\n")
    
  }, error = function(e) {
    cat("Parameter extraction/reconstruction failed:", e$message, "\n")
  })
  
  cat("\n=== Analysis Complete ===\n")
}

# Run the test
test_parameter_formats()