# Debug the matrix dimension error in detail
# Focus on identifying where the "0x1 and 2x1" error originates

library(dirichletprocess)
set.seed(123)

cat("=== DETAILED MATRIX DIMENSION ERROR ANALYSIS ===\n")

# Test with more controlled error reporting
test_initialization_steps <- function(model_name) {
  cat(sprintf("\n--- Debugging %s Model Initialization ---\n", model_name))
  
  # Test data
  x <- matrix(rnorm(20), ncol = 2)
  
  # Step 1: Create mixing distribution
  cat("Step 1: Creating mixing distribution...\n")
  tryCatch({
    md <- MvnormalCreate(list(covModel = model_name))
    cat("   ✓ Mixing distribution created\n")
    
    # Step 2: Create DP object
    cat("Step 2: Creating DP object...\n")
    dp <- DirichletProcessCreate(x, md)
    cat("   ✓ DP object created\n")
    
    # Step 3: Initialize parameters
    cat("Step 3: Initializing parameters...\n")
    
    # Add traceback to see where error occurs
    options(error = function() {
      cat("Error traceback:\n")
      traceback()
    })
    
    dp <- Initialise(dp)
    cat("   ✓ Initialization successful\n")
    
    # Reset error handling
    options(error = NULL)
    
  }, error = function(e) {
    cat(sprintf("   ✗ Error: %s\n", e$message))
    
    # Reset error handling
    options(error = NULL)
  })
}

# Test step-by-step for each model
models_to_test <- c("FULL", "EII", "VII")

for (model in models_to_test) {
  test_initialization_steps(model)
}

# Test with different data sizes
cat("\n--- Testing Different Data Sizes ---\n")
test_data_sizes <- function() {
  model <- "EII"
  
  data_sizes <- list(
    small = matrix(rnorm(6), ncol = 2),     # 3x2
    medium = matrix(rnorm(20), ncol = 2),   # 10x2
    large = matrix(rnorm(100), ncol = 2)    # 50x2
  )
  
  for (size_name in names(data_sizes)) {
    cat(sprintf("Testing %s data (%dx%d)...\n", size_name, 
                nrow(data_sizes[[size_name]]), ncol(data_sizes[[size_name]])))
    
    tryCatch({
      dp <- DirichletProcessCreate(data_sizes[[size_name]], 
                                   MvnormalCreate(list(covModel = model)))
      dp <- Initialise(dp)
      cat("   ✓ Success\n")
      
    }, error = function(e) {
      cat(sprintf("   ✗ Error: %s\n", e$message))
    })
  }
}

test_data_sizes()

# Test with C++ disabled to isolate the issue
cat("\n--- Testing with C++ Completely Disabled ---\n")
test_no_cpp <- function() {
  set_use_cpp(FALSE)
  
  tryCatch({
    x <- matrix(rnorm(20), ncol = 2)
    dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = "EII")))
    dp <- Initialise(dp)
    cat("   ✓ R-only initialization successful\n")
    
    # Try one MCMC iteration
    dp <- Fit(dp, 1, progressBar = FALSE)
    cat("   ✓ R-only MCMC successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ R-only failed: %s\n", e$message))
  })
}

test_no_cpp()

cat("\n=== DETAILED ANALYSIS COMPLETE ===\n")