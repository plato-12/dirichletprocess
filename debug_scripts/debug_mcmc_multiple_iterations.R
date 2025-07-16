# Debug script to investigate MCMC error with multiple iterations

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to investigate MCMC error with multiple iterations
debug_mcmc_multiple_iterations <- function() {
  cat("=== Debugging MCMC Multiple Iterations Error ===\n")
  
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
  
  # Test increasing number of iterations
  for (its in c(1, 2, 5, 10)) {
    cat(sprintf("\n%d. Testing %d MCMC iterations...\n", its, its))
    
    tryCatch({
      dp_obj_test <- Fit(dp_obj, its = its)
      cat(sprintf("  %d MCMC iterations successful.\n", its))
      
    }, error = function(e) {
      cat(sprintf("  %d MCMC iterations failed:", its), e$message, "\n")
      print(e)
      
      # If it fails, let's debug step by step
      if (its > 1) {
        cat(sprintf("    Debugging step-by-step for %d iterations...\n", its))
        
        # Try manual stepping
        dp_temp <- dp_obj
        for (i in 1:its) {
          tryCatch({
            dp_temp <- Fit(dp_temp, its = 1)
            cat(sprintf("      Step %d successful.\n", i))
            
          }, error = function(e2) {
            cat(sprintf("      Step %d failed:", i), e2$message, "\n")
            break
          })
        }
      }
      
      break # Stop testing higher iterations if one fails
    })
  }
  
  cat("\n=== MCMC Multiple Iterations Debug Complete ===\n")
}

# Run the debug
debug_mcmc_multiple_iterations()