# Debug script to trace the exact Lambda + S dimension issue

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to trace Lambda + S dimension issue
debug_lambda_s_dimensions <- function() {
  cat("=== Debugging Lambda + S Dimension Issue ===\n")
  
  # Create simple 1D data - this is the exact data from our test
  test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  
  # Create E model distribution exactly as before
  e_model <- MvnormalCreate(list(
    mu0 = c(0),
    kappa0 = 1,
    nu = 2,
    Lambda = diag(1),
    covModel = "E"
  ))
  
  cat("Test data dimensions:", dim(test_data), "\n")
  cat("E model Lambda dimensions:", dim(e_model$priorParameters$Lambda), "\n")
  cat("E model mu0:", e_model$priorParameters$mu0, "\n")
  
  # Step 1: Call PosteriorParameters directly with test data
  cat("\n1. Testing PosteriorParameters with all test data...\n")
  
  tryCatch({
    post_params <- PosteriorParameters(e_model, test_data)
    cat("PosteriorParameters successful with all data.\n")
    cat("Result t_n dimensions:", dim(post_params$t_n), "\n")
    
  }, error = function(e) {
    cat("PosteriorParameters failed with all data:", e$message, "\n")
    
    # Manual computation to find the issue
    cat("Manual computation:\n")
    
    x <- test_data
    priorParameters <- e_model$priorParameters
    n <- nrow(x)
    d <- ncol(x)
    x_bar <- colMeans(x)
    
    cat("  x dimensions:", dim(x), "\n")
    cat("  n =", n, ", d =", d, "\n")
    cat("  x_bar =", x_bar, ", length =", length(x_bar), "\n")
    cat("  priorParameters$mu0 =", priorParameters$mu0, ", length =", length(priorParameters$mu0), "\n")
    
    # Compute scatter matrix
    if (n > 1) {
      cat("  Computing scatter matrix for n > 1...\n")
      if (e_model$priorParameters$covModel %in% c("E", "V")) {
        cat("  Using E/V model scatter computation...\n")
        S <- (n - 1) * var(x)
        cat("  var(x) =", var(x), ", class =", class(var(x)), "\n")
        S <- matrix(S, 1, 1)
        cat("  S after matrix() =", S, ", dimensions =", dim(S), "\n")
      }
    } else {
      cat("  Using n = 1 scatter matrix...\n")
      S <- matrix(0, d, d)
    }
    
    cat("  Final S dimensions:", dim(S), "\n")
    cat("  priorParameters$Lambda dimensions:", dim(priorParameters$Lambda), "\n")
    
    # Check if Lambda and S are compatible
    tryCatch({
      test_add <- priorParameters$Lambda + S
      cat("  Lambda + S successful, dimensions:", dim(test_add), "\n")
    }, error = function(e2) {
      cat("  Lambda + S failed:", e2$message, "\n")
      cat("  Detailed comparison:\n")
      cat("    Lambda:\n")
      print(priorParameters$Lambda)
      cat("    S:\n")
      print(S)
    })
    
    # Check diff computation
    diff <- x_bar - priorParameters$mu0
    cat("  diff =", diff, ", length =", length(diff), "\n")
    
    # Check outer product
    kappa_n <- priorParameters$kappa0 + n
    outer_term <- (priorParameters$kappa0 * n / kappa_n) * outer(diff, diff)
    cat("  outer_term dimensions:", dim(outer_term), "\n")
  })
  
  cat("\n=== Debug Complete ===\n")
}

# Run the debug
debug_lambda_s_dimensions()