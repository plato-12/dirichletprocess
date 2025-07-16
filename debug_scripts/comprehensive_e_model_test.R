# Comprehensive test suite for E model functionality

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Comprehensive E model test function
comprehensive_e_model_test <- function() {
  cat("=== Comprehensive E Model Test Suite ===\n")
  
  # Test 1: Basic E model creation and initialization
  cat("\n1. Testing E model creation and initialization...\n")
  
  tryCatch({
    # Create 1D data
    test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
    
    # Create E model
    e_model <- MvnormalCreate(list(
      mu0 = c(0),
      kappa0 = 1,
      nu = 2,
      Lambda = diag(1),
      covModel = "E"
    ))
    
    # Test initialization
    dp_e <- DirichletProcessMvnormal(test_data, e_model)
    
    cat("  E model creation and initialization: SUCCESS\n")
    cat("  Initial clusters:", dp_e$numberClusters, "\n")
    cat("  Mu dimensions:", dim(dp_e$clusterParameters$mu), "\n")
    cat("  Sig dimensions:", dim(dp_e$clusterParameters$sig), "\n")
    
  }, error = function(e) {
    cat("  E model creation and initialization: FAILED -", e$message, "\n")
  })
  
  # Test 2: E model MCMC with different data sizes
  cat("\n2. Testing E model MCMC with different data sizes...\n")
  
  for (n in c(5, 10, 20, 50)) {
    cat(sprintf("  Testing with n=%d... ", n))
    
    tryCatch({
      # Generate random 1D data
      test_data <- matrix(rnorm(n), ncol = 1)
      
      e_model <- MvnormalCreate(list(
        mu0 = c(0),
        kappa0 = 1,
        nu = 2,
        Lambda = diag(1),
        covModel = "E"
      ))
      
      dp_e <- DirichletProcessMvnormal(test_data, e_model)
      dp_e_mcmc <- Fit(dp_e, its = 10)
      
      cat("SUCCESS\n")
      
    }, error = function(e) {
      cat("FAILED -", e$message, "\n")
    })
  }
  
  # Test 3: E model with different prior parameters
  cat("\n3. Testing E model with different prior parameters...\n")
  
  priors_list <- list(
    list(mu0 = c(0), kappa0 = 1, nu = 2, Lambda = diag(1)),
    list(mu0 = c(5), kappa0 = 0.1, nu = 3, Lambda = diag(1) * 2),
    list(mu0 = c(-2), kappa0 = 10, nu = 5, Lambda = diag(1) * 0.5)
  )
  
  for (i in seq_along(priors_list)) {
    cat(sprintf("  Testing prior set %d... ", i))
    
    tryCatch({
      test_data <- matrix(rnorm(10), ncol = 1)
      
      priors <- priors_list[[i]]
      priors$covModel <- "E"
      
      e_model <- MvnormalCreate(priors)
      dp_e <- DirichletProcessMvnormal(test_data, e_model)
      dp_e_mcmc <- Fit(dp_e, its = 10)
      
      cat("SUCCESS\n")
      
    }, error = function(e) {
      cat("FAILED -", e$message, "\n")
    })
  }
  
  # Test 4: E model with different initial cluster numbers
  cat("\n4. Testing E model with different initial cluster numbers...\n")
  
  for (k in c(1, 2, 3, 5)) {
    cat(sprintf("  Testing with %d initial clusters... ", k))
    
    tryCatch({
      test_data <- matrix(rnorm(20), ncol = 1)
      
      e_model <- MvnormalCreate(list(
        mu0 = c(0),
        kappa0 = 1,
        nu = 2,
        Lambda = diag(1),
        covModel = "E"
      ))
      
      dp_e <- DirichletProcessMvnormal(test_data, e_model)
      dp_e <- Initialise(dp_e, numInitialClusters = k)
      dp_e_mcmc <- Fit(dp_e, its = 10)
      
      cat("SUCCESS\n")
      
    }, error = function(e) {
      cat("FAILED -", e$message, "\n")
    })
  }
  
  # Test 5: E model clustering correctness
  cat("\n5. Testing E model clustering correctness...\n")
  
  tryCatch({
    # Create data with clear clusters
    set.seed(123)
    data1 <- rnorm(10, mean = 0, sd = 0.5)
    data2 <- rnorm(10, mean = 5, sd = 0.5)
    test_data <- matrix(c(data1, data2), ncol = 1)
    
    e_model <- MvnormalCreate(list(
      mu0 = c(2.5),
      kappa0 = 0.1,
      nu = 3,
      Lambda = diag(1),
      covModel = "E"
    ))
    
    dp_e <- DirichletProcessMvnormal(test_data, e_model)
    dp_e_mcmc <- Fit(dp_e, its = 100)
    
    # Check if it finds reasonable clusters
    unique_clusters <- length(unique(dp_e_mcmc$clusterLabels))
    cat(sprintf("  Found %d unique clusters (expected around 2): ", unique_clusters))
    
    if (unique_clusters >= 2 && unique_clusters <= 5) {
      cat("SUCCESS\n")
    } else {
      cat("QUESTIONABLE\n")
    }
    
  }, error = function(e) {
    cat("  Clustering correctness: FAILED -", e$message, "\n")
  })
  
  # Test 6: E model parameter estimation
  cat("\n6. Testing E model parameter estimation...\n")
  
  tryCatch({
    # Create data with known parameters
    set.seed(456)
    true_mean <- 3.0
    true_var <- 2.0
    test_data <- matrix(rnorm(50, mean = true_mean, sd = sqrt(true_var)), ncol = 1)
    
    e_model <- MvnormalCreate(list(
      mu0 = c(0),
      kappa0 = 0.1,
      nu = 3,
      Lambda = diag(1),
      covModel = "E"
    ))
    
    dp_e <- DirichletProcessMvnormal(test_data, e_model)
    dp_e_mcmc <- Fit(dp_e, its = 100)
    
    # Get parameter estimates (assuming single cluster)
    if (dp_e_mcmc$numberClusters == 1) {
      est_mean <- dp_e_mcmc$clusterParameters$mu[1, 1, 1]
      est_var <- dp_e_mcmc$clusterParameters$sig[1, 1]
      
      cat(sprintf("  True mean: %.2f, Estimated: %.2f\n", true_mean, est_mean))
      cat(sprintf("  True var: %.2f, Estimated: %.2f\n", true_var, est_var))
      
      if (abs(est_mean - true_mean) < 1.0 && abs(est_var - true_var) < 3.0) {
        cat("  Parameter estimation: SUCCESS\n")
      } else {
        cat("  Parameter estimation: QUESTIONABLE\n")
      }
    } else {
      cat("  Parameter estimation: SKIPPED (multiple clusters)\n")
    }
    
  }, error = function(e) {
    cat("  Parameter estimation: FAILED -", e$message, "\n")
  })
  
  cat("\n=== Comprehensive E Model Test Complete ===\n")
}

# Run the comprehensive test
comprehensive_e_model_test()