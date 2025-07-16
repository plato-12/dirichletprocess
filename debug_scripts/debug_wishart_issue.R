# Debug the rWishart numerical stability issue
# Focus on the BLAS/LAPACK error code -4

library(dirichletprocess)
set.seed(123)

cat("=== WISHART NUMERICAL STABILITY ANALYSIS ===\n")

# Test rWishart with different parameters
test_wishart_parameters <- function() {
  cat("Testing rWishart with different parameter combinations...\n")
  
  # Test simple case first
  cat("1. Testing simple 2x2 identity matrix...\n")
  tryCatch({
    result <- rWishart(1, 3, diag(2))
    cat("   ✓ Simple case successful\n")
  }, error = function(e) {
    cat(sprintf("   ✗ Simple case failed: %s\n", e$message))
  })
  
  # Test with different degrees of freedom
  cat("2. Testing different degrees of freedom...\n")
  for (nu in c(2, 3, 4, 5)) {
    tryCatch({
      result <- rWishart(1, nu, diag(2))
      cat(sprintf("   ✓ nu=%d successful\n", nu))
    }, error = function(e) {
      cat(sprintf("   ✗ nu=%d failed: %s\n", nu, e$message))
    })
  }
  
  # Test with different Lambda matrices
  cat("3. Testing different Lambda matrices...\n")
  
  # Test cases
  lambda_cases <- list(
    identity = diag(2),
    scaled = diag(2) * 0.1,
    small = diag(2) * 0.01,
    large = diag(2) * 10
  )
  
  for (case_name in names(lambda_cases)) {
    Lambda <- lambda_cases[[case_name]]
    
    cat(sprintf("   Testing %s matrix...\n", case_name))
    
    # Check matrix properties
    det_Lambda <- det(Lambda)
    eigenvals <- eigen(Lambda)$values
    
    cat(sprintf("     Det(Lambda): %g\n", det_Lambda))
    cat(sprintf("     Eigenvalues: %s\n", paste(round(eigenvals, 6), collapse = ", ")))
    cat(sprintf("     Positive definite: %s\n", all(eigenvals > 0)))
    
    tryCatch({
      result <- rWishart(1, 3, Lambda)
      cat(sprintf("     ✓ %s successful\n", case_name))
    }, error = function(e) {
      cat(sprintf("     ✗ %s failed: %s\n", case_name, e$message))
    })
  }
}

test_wishart_parameters()

# Test the actual prior parameters from MVNormal
cat("\n--- Testing Actual MVNormal Prior Parameters ---\n")
test_mvnormal_priors <- function() {
  models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
  
  for (model in models) {
    cat(sprintf("Testing %s model prior parameters...\n", model))
    
    tryCatch({
      # Create mixing distribution
      md <- MvnormalCreate(list(covModel = model))
      prior_params <- md$priorParameters
      
      cat(sprintf("   nu: %g\n", prior_params$nu))
      cat(sprintf("   Lambda det: %g\n", det(prior_params$Lambda)))
      
      # Test rWishart with these parameters
      result <- rWishart(1, prior_params$nu, prior_params$Lambda)
      cat(sprintf("   ✓ %s prior rWishart successful\n", model))
      
    }, error = function(e) {
      cat(sprintf("   ✗ %s prior rWishart failed: %s\n", model, e$message))
    })
  }
}

test_mvnormal_priors()

# Test with actual data posterior parameters
cat("\n--- Testing Posterior Parameters ---\n")
test_posterior_parameters <- function() {
  x <- matrix(rnorm(20), ncol = 2)
  
  tryCatch({
    # Create a simple mixing distribution
    md <- MvnormalCreate(list(covModel = "EII"))
    
    # Compute posterior parameters manually
    post_params <- PosteriorParameters(md, x)
    
    cat("Posterior parameters computed:\n")
    cat(sprintf("   nu_n: %g\n", post_params$nu_n))
    cat(sprintf("   t_n det: %g\n", det(post_params$t_n)))
    cat(sprintf("   t_n eigenvalues: %s\n", 
                paste(round(eigen(post_params$t_n)$values, 6), collapse = ", ")))
    
    # Test rWishart with posterior parameters
    result <- rWishart(1, post_params$nu_n, post_params$t_n)
    cat("   ✓ Posterior rWishart successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ Posterior rWishart failed: %s\n", e$message))
  })
}

test_posterior_parameters()

# Test matrix regularization
cat("\n--- Testing Matrix Regularization ---\n")
test_regularization <- function() {
  x <- matrix(rnorm(20), ncol = 2)
  
  tryCatch({
    md <- MvnormalCreate(list(covModel = "EII"))
    post_params <- PosteriorParameters(md, x)
    
    # Add regularization to t_n
    t_n_reg <- post_params$t_n + diag(nrow(post_params$t_n)) * 1e-6
    
    cat("Regularized matrix:\n")
    cat(sprintf("   Det: %g\n", det(t_n_reg)))
    cat(sprintf("   Eigenvalues: %s\n", 
                paste(round(eigen(t_n_reg)$values, 6), collapse = ", ")))
    
    # Test rWishart with regularized matrix
    result <- rWishart(1, post_params$nu_n, t_n_reg)
    cat("   ✓ Regularized rWishart successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ Regularized rWishart failed: %s\n", e$message))
  })
}

test_regularization()

cat("\n=== WISHART ANALYSIS COMPLETE ===\n")