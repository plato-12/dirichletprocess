# Isolate the specific atime integration issue

library(dirichletprocess)
devtools::load_all()

# Load atime if available
if (!require(atime, quietly = TRUE)) {
  cat("atime package not available\n")
  quit()
}

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Isolate atime issue
isolate_atime_issue <- function() {
  cat("=== Isolating atime Integration Issue ===\n")
  
  # Test the exact error pattern mentioned: "values must be length 2"
  cat("\n1. Testing simple case that works...\n")
  
  tryCatch({
    # This should work - returns consistent format
    simple_test <- atime::atime(
      N = c(5, 10),
      setup = {
        data <- rnorm(N)
      },
      test = {
        # Always return same length vector
        c(mean(data), sd(data))
      }
    )
    cat("  Simple consistent return test: SUCCESS\n")
    
  }, error = function(e) {
    cat("  Simple test failed:", e$message, "\n")
  })
  
  # Test 2: Inconsistent return lengths (this should fail)
  cat("\n2. Testing inconsistent return lengths...\n")
  
  tryCatch({
    # This might fail - inconsistent return format
    inconsistent_test <- atime::atime(
      N = c(5, 10),
      setup = {
        data <- rnorm(N)
      },
      test = {
        # Return different lengths based on data
        if (N < 8) {
          mean(data)  # Returns length 1
        } else {
          c(mean(data), sd(data))  # Returns length 2
        }
      }
    )
    cat("  Inconsistent return test: UNEXPECTED SUCCESS\n")
    
  }, error = function(e) {
    cat("  Inconsistent return test failed (expected):", e$message, "\n")
  })
  
  # Test 3: Test with dirichletprocess objects
  cat("\n3. Testing DP object return formats...\n")
  
  tryCatch({
    # Test DP with consistent return format
    dp_test <- atime::atime(
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
      test = {
        md <- MvnormalCreate(prior_params)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        
        # Return consistent format
        c(dp$numberClusters, length(dp$clusterLabels))
      }
    )
    cat("  DP consistent return test: SUCCESS\n")
    
  }, error = function(e) {
    cat("  DP test failed:", e$message, "\n")
    print(e)
  })
  
  # Test 4: Test potential issue with DP object returns
  cat("\n4. Testing DP object direct return (potential issue)...\n")
  
  tryCatch({
    # This might be the issue - returning complex objects
    dp_object_test <- atime::atime(
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
      test = {
        md <- MvnormalCreate(prior_params)
        dp <- DirichletProcessCreate(test_data, md)
        dp <- Initialise(dp, numInitialClusters = 1)
        
        # Return the DP object directly (might cause issues)
        dp
      }
    )
    cat("  DP object return test: SUCCESS\n")
    
  }, error = function(e) {
    cat("  DP object return test failed:", e$message, "\n")
    print(e)
  })
  
  cat("\n=== Isolation Test Complete ===\n")
}

# Run the test
isolate_atime_issue()