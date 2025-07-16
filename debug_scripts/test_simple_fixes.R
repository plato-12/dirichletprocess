# Test simple fixes
library(dirichletprocess)
set_use_cpp(FALSE)  # Use R implementation to test safe rWishart

# Source the safe rWishart file
source("R/safe_wishart.R")

cat("=== TESTING SIMPLE FIXES ===\n")

# Test 1: Safe rWishart wrapper
cat("Testing safe rWishart wrapper...\n")
tryCatch({
  # Test with well-conditioned matrix
  Lambda <- diag(2)
  result <- safe_rWishart(1, 3, Lambda)
  cat("  ✓ Well-conditioned matrix: SUCCESS\n")
  
  # Test with slightly ill-conditioned matrix
  Lambda_ill <- matrix(c(1, 0.999, 0.999, 1), 2, 2)
  result <- safe_rWishart(1, 3, Lambda_ill)
  cat("  ✓ Ill-conditioned matrix: SUCCESS\n")
  
  # Test matrix analysis
  analysis <- analyze_matrix_conditioning(Lambda_ill)
  cat(sprintf("  Matrix condition number: %.2e\n", analysis$condition_number))
  cat(sprintf("  Needs regularization: %s\n", analysis$needs_regularization))
  
}, error = function(e) {
  cat(sprintf("  ✗ Safe rWishart failed: %s\n", e$message))
})

# Test 2: Constrained models work
cat("\nTesting constrained models...\n")
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
test_data <- matrix(rnorm(40), ncol = 2)

for (model in models) {
  tryCatch({
    # Create model
    md <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = model
    ))
    
    # Test basic functionality
    dp <- DirichletProcessCreate(test_data, md)
    dp <- Initialise(dp, numInitialClusters = 2)
    dp <- Fit(dp, 10, progressBar = FALSE)
    
    cat(sprintf("  ✓ %s: SUCCESS (%d clusters)\n", model, dp$numberClusters))
  }, error = function(e) {
    cat(sprintf("  ✗ %s: FAILED - %s\n", model, e$message))
  })
}

cat("\n=== SIMPLE FIXES TEST COMPLETE ===\n")