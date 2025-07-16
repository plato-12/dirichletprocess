# Debug script to investigate MCMC "incorrect number of subscripts" error

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to investigate MCMC error
debug_mcmc_subscripts <- function() {
  cat("=== Debugging MCMC Subscripts Error ===\n")
  
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
  
  # Initialize the DP object
  dp_obj <- DirichletProcessMvnormal(test_data, e_model)
  
  cat("DP object created successfully.\n")
  cat("Initial cluster parameters:\n")
  cat("  mu dimensions:", dim(dp_obj$clusterParameters$mu), "\n")
  cat("  sig dimensions:", dim(dp_obj$clusterParameters$sig), "\n")
  
  # Test single MCMC step
  cat("\n1. Testing single MCMC iteration...\n")
  
  tryCatch({
    # Try to run a single iteration
    dp_obj_step <- Fit(dp_obj, its = 1)
    cat("Single MCMC iteration successful.\n")
    
  }, error = function(e) {
    cat("Single MCMC iteration failed:", e$message, "\n")
    print(e)
    
    # Let's test the individual MCMC components
    cat("\n  Testing individual MCMC components...\n")
    
    # Test ClusterComponentUpdate
    cat("  Testing ClusterComponentUpdate...\n")
    tryCatch({
      dp_obj_cluster <- ClusterComponentUpdate(dp_obj)
      cat("  ClusterComponentUpdate successful.\n")
      
    }, error = function(e2) {
      cat("  ClusterComponentUpdate failed:", e2$message, "\n")
    })
    
    # Test ClusterParameterUpdate
    cat("  Testing ClusterParameterUpdate...\n")
    tryCatch({
      dp_obj_params <- ClusterParameterUpdate(dp_obj)
      cat("  ClusterParameterUpdate successful.\n")
      
    }, error = function(e3) {
      cat("  ClusterParameterUpdate failed:", e3$message, "\n")
      print(e3)
    })
    
    # Test UpdateAlpha
    cat("  Testing UpdateAlpha...\n")
    tryCatch({
      dp_obj_alpha <- UpdateAlpha(dp_obj)
      cat("  UpdateAlpha successful.\n")
      
    }, error = function(e4) {
      cat("  UpdateAlpha failed:", e4$message, "\n")
    })
  })
  
  # Test with verbose output to see where it fails
  cat("\n2. Testing with verbose output...\n")
  
  tryCatch({
    dp_obj_verbose <- Fit(dp_obj, its = 1, verbose = TRUE)
    cat("Verbose MCMC successful.\n")
    
  }, error = function(e) {
    cat("Verbose MCMC failed:", e$message, "\n")
  })
  
  cat("\n=== MCMC Debug Complete ===\n")
}

# Run the debug
debug_mcmc_subscripts()