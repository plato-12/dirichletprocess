# Debug script to analyze PosteriorParameters issue with E/V models
# This script traces the exact dimension mismatch

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to trace PosteriorParameters
debug_posterior_parameters <- function() {
  cat("=== Debugging PosteriorParameters for E/V Models ===\n")
  
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
  
  cat("E model created successfully.\n")
  cat("E model priorParameters:\n")
  str(e_model$priorParameters)
  
  # Test with different data sizes
  cat("\n1. Testing with single observation...\n")
  
  single_data <- matrix(c(1), ncol = 1)
  
  tryCatch({
    post_params <- PosteriorParameters(e_model, single_data)
    cat("PosteriorParameters successful with single observation.\n")
    
  }, error = function(e) {
    cat("PosteriorParameters failed with single observation:", e$message, "\n")
    
    # Let's manually trace what happens
    cat("Manual trace:\n")
    
    x <- single_data
    priorParameters <- e_model$priorParameters
    n <- nrow(x)
    d <- ncol(x)
    x_bar <- colMeans(x)
    
    cat("  n =", n, "\n")
    cat("  d =", d, "\n")
    cat("  x_bar =", x_bar, "\n")
    
    # Compute scatter matrix
    if (n > 1) {
      cat("  n > 1, computing scatter matrix\n")
      if (e_model$priorParameters$covModel %in% c("E", "V")) {
        S <- (n - 1) * var(x)
        S <- matrix(S, 1, 1)
      }
    } else {
      cat("  n = 1, using zero scatter matrix\n")
      S <- matrix(0, d, d)
    }
    
    cat("  S dimensions:", dim(S), "\n")
    cat("  S:\n")
    print(S)
    
    cat("  priorParameters$Lambda dimensions:", dim(priorParameters$Lambda), "\n")
    cat("  priorParameters$Lambda:\n")
    print(priorParameters$Lambda)
    
    # Check if we can add them
    tryCatch({
      test_add <- priorParameters$Lambda + S
      cat("  Lambda + S successful, dimensions:", dim(test_add), "\n")
    }, error = function(e2) {
      cat("  Lambda + S failed:", e2$message, "\n")
    })
    
    # Compute the diff term
    diff <- x_bar - priorParameters$mu0
    cat("  diff =", diff, "\n")
    
    # Compute outer product
    kappa_n <- priorParameters$kappa0 + n
    outer_term <- (priorParameters$kappa0 * n / kappa_n) * outer(diff, diff)
    cat("  outer_term dimensions:", dim(outer_term), "\n")
    cat("  outer_term:\n")
    print(outer_term)
    
    # Try the full computation
    tryCatch({
      t_n <- priorParameters$Lambda + S + outer_term
      cat("  Full computation successful, t_n dimensions:", dim(t_n), "\n")
    }, error = function(e3) {
      cat("  Full computation failed:", e3$message, "\n")
    })
  })
  
  cat("\n2. Testing with multiple observations...\n")
  
  multi_data <- matrix(c(1, 2, 3), ncol = 1)
  
  tryCatch({
    post_params <- PosteriorParameters(e_model, multi_data)
    cat("PosteriorParameters successful with multiple observations.\n")
    
  }, error = function(e) {
    cat("PosteriorParameters failed with multiple observations:", e$message, "\n")
  })
  
  cat("\n=== Debug Complete ===\n")
}

# Run the debug
debug_posterior_parameters()