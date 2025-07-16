# Diagnostic script to isolate C++ vs R behavior for constrained models
# This addresses Issue 1: C++ Implementation Matrix Dimension Bugs

library(dirichletprocess)
set.seed(123)

# Create test data
x <- matrix(rnorm(20), ncol = 2)

cat("=== DIAGNOSTIC: C++ vs R Behavior Analysis ===\n")

# Test function to isolate C++ vs R dispatch
test_cpp_vs_r <- function(model_name) {
  cat(sprintf("\n--- Testing %s Model ---\n", model_name))
  
  # Test 1: C++ enabled
  cat("1. Testing with C++ enabled:\n")
  enable_cpp_samplers()
  set_use_cpp(TRUE)
  
  tryCatch({
    dp_cpp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model_name)))
    dp_cpp <- Initialise(dp_cpp)
    cat("   ✓ C++ initialization successful\n")
    
    # Try one MCMC iteration
    dp_cpp <- Fit(dp_cpp, 1, progressBar = FALSE)
    cat("   ✓ C++ MCMC iteration successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ C++ failed: %s\n", e$message))
  })
  
  # Test 2: R fallback only
  cat("2. Testing with R fallback only:\n")
  set_use_cpp(FALSE)
  
  tryCatch({
    dp_r <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model_name)))
    dp_r <- Initialise(dp_r)
    cat("   ✓ R initialization successful\n")
    
    # Try one MCMC iteration
    dp_r <- Fit(dp_r, 1, progressBar = FALSE)
    cat("   ✓ R MCMC iteration successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ R failed: %s\n", e$message))
  })
  
  # Test 3: Individual C++ functions
  cat("3. Testing individual C++ functions:\n")
  enable_cpp_samplers()
  set_use_cpp(TRUE)
  
  tryCatch({
    # Test if C++ functions exist
    pkg_ns <- getNamespace("dirichletprocess")
    cpp_functions <- c(
      "mvnormal_posterior_draw_cpp",
      "mvnormal_posterior_parameters_cpp",
      "mvnormal_prior_draw_cpp",
      "mvnormal_likelihood_cpp"
    )
    
    for (func_name in cpp_functions) {
      if (exists(func_name, where = pkg_ns)) {
        cat(sprintf("   ✓ %s exists\n", func_name))
      } else {
        cat(sprintf("   ✗ %s missing\n", func_name))
      }
    }
    
  }, error = function(e) {
    cat(sprintf("   ✗ C++ function check failed: %s\n", e$message))
  })
}

# Test all constrained models
constrained_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")

for (model in constrained_models) {
  test_cpp_vs_r(model)
}

# Test FULL model for comparison
cat("\n--- Testing FULL Model (Reference) ---\n")
test_cpp_vs_r("FULL")

cat("\n=== DIAGNOSTIC COMPLETE ===\n")