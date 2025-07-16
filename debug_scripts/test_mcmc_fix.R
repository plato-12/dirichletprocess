# Test script to verify the fix for MCMC subscripts error

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Test function to verify the MCMC fix
test_mcmc_fix <- function() {
  cat("=== Testing MCMC Fix for E/V Models ===\n")
  
  # Create simple 1D data
  test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  
  # Test E model
  cat("\n1. Testing E model MCMC...\n")
  
  tryCatch({
    e_model <- MvnormalCreate(list(
      mu0 = c(0),
      kappa0 = 1,
      nu = 2,
      Lambda = diag(1),
      covModel = "E"
    ))
    
    dp_e <- DirichletProcessMvnormal(test_data, e_model)
    
    # Test with increasing iterations
    for (its in c(1, 2, 5, 10, 20)) {
      cat(sprintf("  Testing %d iterations... ", its))
      
      dp_e_mcmc <- Fit(dp_e, its = its)
      cat("SUCCESS\n")
    }
    
    cat("E model MCMC fully successful!\n")
    
  }, error = function(e) {
    cat("E model MCMC failed:", e$message, "\n")
  })
  
  # Test V model
  cat("\n2. Testing V model MCMC...\n")
  
  tryCatch({
    v_model <- MvnormalCreate(list(
      mu0 = c(0),
      kappa0 = 1,
      nu = 2,
      Lambda = diag(1),
      covModel = "V"
    ))
    
    dp_v <- DirichletProcessMvnormal(test_data, v_model)
    
    # Test with increasing iterations
    for (its in c(1, 2, 5, 10, 20)) {
      cat(sprintf("  Testing %d iterations... ", its))
      
      dp_v_mcmc <- Fit(dp_v, its = its)
      cat("SUCCESS\n")
    }
    
    cat("V model MCMC fully successful!\n")
    
  }, error = function(e) {
    cat("V model MCMC failed:", e$message, "\n")
  })
  
  # Test multivariate models to ensure no regression
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
    
    dp_eii <- DirichletProcessMvnormal(test_data_2d, eii_model)
    dp_eii_mcmc <- Fit(dp_eii, its = 10)
    
    cat("EII model MCMC successful!\n")
    
    # FULL model
    full_model <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = "FULL"
    ))
    
    dp_full <- DirichletProcessMvnormal(test_data_2d, full_model)
    dp_full_mcmc <- Fit(dp_full, its = 10)
    
    cat("FULL model MCMC successful!\n")
    
  }, error = function(e) {
    cat("Multivariate models failed:", e$message, "\n")
  })
  
  cat("\n=== MCMC Fix Test Complete ===\n")
}

# Run the test
test_mcmc_fix()