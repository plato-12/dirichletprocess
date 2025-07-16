# Minimal atime benchmark to identify the exact issue

library(dirichletprocess)
devtools::load_all()

# Load atime if available
if (!require(atime, quietly = TRUE)) {
  cat("atime package not available\n")
  quit()
}

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Minimal benchmark test
minimal_atime_benchmark <- function() {
  cat("=== Minimal atime Benchmark Test ===\n")
  
  # Create helper functions exactly like in the main benchmark
  create_prior_parameters <- function(dimensions, model_name) {
    # Handle univariate models
    if (model_name %in% c("E", "V") && dimensions > 1) {
      stop("Models E and V are only for univariate data (dimensions = 1)")
    }
    
    # Base parameters
    base_params <- list(
      mu0 = rep(0, dimensions),
      kappa0 = 1,
      nu = dimensions + 1,
      Lambda = diag(dimensions),
      covModel = model_name
    )
    
    return(base_params)
  }
  
  # Test 1: Simple benchmark with minimal settings
  cat("\n1. Testing minimal benchmark setup...\n")
  
  tryCatch({
    # Very simple benchmark - fewer iterations, smaller data
    simple_result <- atime::atime(
      N = 2^c(3, 4),  # Just 8 and 16 samples
      setup = {
        # Simple data
        test_data <- matrix(rnorm(N * 2), ncol = 2)
        prior_FULL <- create_prior_parameters(2, "FULL")
        prior_EII <- create_prior_parameters(2, "EII")
      },
      
      FULL = {
        md <- MvnormalCreate(prior_FULL)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        result <- Fit(dp, 2, progressBar = FALSE)  # Only 2 MCMC iterations
        # Return simple metrics instead of complex object
        result$numberClusters
      },
      
      EII = {
        md <- MvnormalCreate(prior_EII)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        result <- Fit(dp, 2, progressBar = FALSE)  # Only 2 MCMC iterations
        # Return simple metrics instead of complex object
        result$numberClusters
      },
      
      times = 2  # Only 2 repetitions
    )
    
    cat("  Minimal benchmark: SUCCESS\n")
    cat("  Result structure:\n")
    str(simple_result, max.level = 2)
    
  }, error = function(e) {
    cat("  Minimal benchmark failed:", e$message, "\n")
    print(e)
  })
  
  # Test 2: Test return value issues specifically
  cat("\n2. Testing return value consistency...\n")
  
  tryCatch({
    # Test different return value formats
    return_test <- atime::atime(
      N = c(5, 10),
      setup = {
        test_data <- matrix(rnorm(N * 2), ncol = 2)
        prior_params <- create_prior_parameters(2, "FULL")
      },
      
      single_value = {
        md <- MvnormalCreate(prior_params)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        # Return single value
        dp$numberClusters
      },
      
      vector_value = {
        md <- MvnormalCreate(prior_params)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        # Return consistent vector
        c(dp$numberClusters, length(dp$clusterLabels))
      },
      
      times = 2
    )
    
    cat("  Return value test: SUCCESS\n")
    
  }, error = function(e) {
    cat("  Return value test failed:", e$message, "\n")
    print(e)
  })
  
  # Test 3: Test the specific error pattern mentioned
  cat("\n3. Testing problematic patterns...\n")
  
  tryCatch({
    # This might trigger the "length 2 vs length 1" error
    problematic_test <- atime::atime(
      N = c(5, 10),
      setup = {
        test_data <- matrix(rnorm(N * 2), ncol = 2)
        prior_params <- create_prior_parameters(2, "FULL")
      },
      
      # This might cause issues if it returns different structures
      variable_return = {
        md <- MvnormalCreate(prior_params)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        result <- Fit(dp, 1, progressBar = FALSE)
        
        # Return the whole DP object (might cause format issues)
        result
      },
      
      times = 2
    )
    
    cat("  Problematic pattern test: UNEXPECTED SUCCESS\n")
    
  }, error = function(e) {
    cat("  Problematic pattern test failed (this might be the issue):", e$message, "\n")
    print(e)
  })
  
  cat("\n=== Minimal Benchmark Test Complete ===\n")
}

# Run the test
minimal_atime_benchmark()