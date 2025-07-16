# Test a simple fix by providing explicit parameters
library(dirichletprocess)
set.seed(123)

cat("=== TESTING SIMPLE FIX ===\n")

# Test the working case
test_explicit_parameters <- function() {
  cat("Testing with explicit parameters...\n")
  
  x <- matrix(rnorm(20), ncol = 2)
  
  tryCatch({
    # Create with explicit parameters
    md <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = "EII"
    ))
    
    dp <- DirichletProcessCreate(x, md)
    dp <- Initialise(dp)
    
    cat("   ✓ Explicit parameters work\n")
    
    # Try a short MCMC run
    dp <- Fit(dp, 5, progressBar = FALSE)
    cat("   ✓ Short MCMC run successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ Explicit parameters failed: %s\n", e$message))
  })
}

test_explicit_parameters()

# Test all constrained models with explicit parameters
cat("\n--- Testing All Constrained Models with Explicit Parameters ---\n")
test_all_constrained <- function() {
  models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  x <- matrix(rnorm(20), ncol = 2)
  
  for (model in models) {
    cat(sprintf("Testing %s...\n", model))
    
    tryCatch({
      md <- MvnormalCreate(list(
        mu0 = c(0, 0),
        kappa0 = 1,
        nu = 3,
        Lambda = diag(2),
        covModel = model
      ))
      
      dp <- DirichletProcessCreate(x, md)
      dp <- Initialise(dp)
      
      # Try a short MCMC run
      dp <- Fit(dp, 5, progressBar = FALSE)
      cat(sprintf("   ✓ %s successful\n", model))
      
    }, error = function(e) {
      cat(sprintf("   ✗ %s failed: %s\n", model, e$message))
    })
  }
}

test_all_constrained()

# Test the benchmark validation function
cat("\n--- Testing Benchmark Validation ---\n")
test_benchmark_validation <- function() {
  tryCatch({
    source("benchmark/atime/benchmark-covariance-models-comprehensive.R")
    
    # Run the validation
    result <- validate_all_models()
    cat(sprintf("   Validation result: %s\n", result))
    
  }, error = function(e) {
    cat(sprintf("   ✗ Benchmark validation failed: %s\n", e$message))
  })
}

test_benchmark_validation()

cat("\n=== SIMPLE FIX TEST COMPLETE ===\n")