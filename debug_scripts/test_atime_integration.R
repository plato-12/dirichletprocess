# Test script to analyze atime benchmark integration issues

library(dirichletprocess)
devtools::load_all()

# Load atime if available
if (!require(atime, quietly = TRUE)) {
  cat("atime package not available, installing...\n")
  install.packages("atime")
  library(atime)
}

# Test atime integration issues
test_atime_integration <- function() {
  cat("=== Analyzing atime Framework Integration Issues ===\n")
  
  # Test 1: Simple atime test to understand expected format
  cat("\n1. Testing basic atime functionality...\n")
  
  tryCatch({
    # Simple test function that returns different types
    test_result <- atime::atime(
      N = c(10, 20),
      setup = {
        test_data <- rnorm(N)
      },
      single_value = {
        mean(test_data)
      },
      multiple_values = {
        c(mean(test_data), sd(test_data))
      }
    )
    
    cat("  Basic atime test successful.\n")
    cat("  Result structure:\n")
    str(test_result)
    
  }, error = function(e) {
    cat("  Basic atime test failed:", e$message, "\n")
  })
  
  # Test 2: Test with dirichletprocess objects
  cat("\n2. Testing dirichletprocess with atime...\n")
  
  tryCatch({
    # Test if the issue is with return value format
    dp_result <- atime::atime(
      N = c(5, 10),
      setup = {
        test_data <- matrix(rnorm(N * 2), ncol = 2)
        prior_params <- list(
          mu0 = c(0, 0),
          kappa0 = 1,
          nu = 3,
          Lambda = diag(2),
          covModel = "FULL"
        )
      },
      dp_create = {
        md <- MvnormalCreate(prior_params)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        # Return a simple metric instead of the whole object
        dp$numberClusters
      }
    )
    
    cat("  Dirichletprocess atime test successful.\n")
    cat("  Result structure:\n")
    str(dp_result)
    
  }, error = function(e) {
    cat("  Dirichletprocess atime test failed:", e$message, "\n")
    print(e)
  })
  
  # Test 3: Test return value issues mentioned in action plan
  cat("\n3. Testing return value format issues...\n")
  
  tryCatch({
    # Test the exact error pattern: "values must be length 2, but FUN(X[[1]]) result is length 1"
    problematic_result <- atime::atime(
      N = c(5, 10),
      setup = {
        test_data <- matrix(rnorm(N), ncol = 1)
      },
      single_return = {
        # This might cause the "length 1" issue
        length(test_data)
      },
      double_return = {
        # This should work fine
        c(length(test_data), mean(test_data))
      }
    )
    
    cat("  Return value test successful.\n")
    
  }, error = function(e) {
    cat("  Return value test failed:", e$message, "\n")
    print(e)
  })
  
  # Test 4: Test with benchmark helper functions
  cat("\n4. Testing benchmark helper functions...\n")
  
  # Check if helper functions exist
  source_file <- "benchmark/atime/benchmark-covariance-models-comprehensive.R"
  if (file.exists(source_file)) {
    tryCatch({
      source(source_file)
      
      # Test create_prior_parameters function
      prior_test <- create_prior_parameters(2, "FULL")
      cat("  create_prior_parameters works.\n")
      
      # Test run_atime_benchmark function
      cat("  Testing run_atime_benchmark...\n")
      benchmark_result <- run_atime_benchmark()
      cat("  run_atime_benchmark successful.\n")
      
    }, error = function(e) {
      cat("  Benchmark functions failed:", e$message, "\n")
      print(e)
    })
  } else {
    cat("  Benchmark file not found at:", source_file, "\n")
  }
  
  cat("\n=== atime Integration Analysis Complete ===\n")
}

# Run the test
test_atime_integration()