# Debug parameter creation in MvnormalCreate
# Focus on why Lambda is NULL or has wrong dimensions

library(dirichletprocess)
set.seed(123)

cat("=== PARAMETER CREATION DEBUG ===\n")

# Test MvnormalCreate parameter creation
debug_mvnormal_create <- function(model_name) {
  cat(sprintf("\n--- Debugging %s Parameter Creation ---\n", model_name))
  
  tryCatch({
    # Step-by-step parameter creation
    cat("1. Creating with minimal parameters...\n")
    
    # Test different parameter specifications
    param_specs <- list(
      minimal = list(covModel = model_name),
      explicit = list(
        mu0 = c(0, 0),
        kappa0 = 1,
        nu = 3,
        Lambda = diag(2),
        covModel = model_name
      )
    )
    
    for (spec_name in names(param_specs)) {
      cat(sprintf("   Testing %s specification...\n", spec_name))
      
      md <- MvnormalCreate(param_specs[[spec_name]])
      
      # Check all parameters
      cat(sprintf("     mu0: %s\n", paste(md$priorParameters$mu0, collapse = ", ")))
      cat(sprintf("     kappa0: %g\n", md$priorParameters$kappa0))
      cat(sprintf("     nu: %g\n", md$priorParameters$nu))
      cat(sprintf("     Lambda class: %s\n", class(md$priorParameters$Lambda)))
      cat(sprintf("     Lambda dim: %s\n", 
                  if (is.null(md$priorParameters$Lambda)) "NULL" 
                  else paste(dim(md$priorParameters$Lambda), collapse = "x")))
      cat(sprintf("     covModel: %s\n", md$priorParameters$covModel))
      
      # Test if Lambda is usable
      if (!is.null(md$priorParameters$Lambda)) {
        cat(sprintf("     Lambda det: %g\n", det(md$priorParameters$Lambda)))
        cat(sprintf("     Lambda eigenvalues: %s\n", 
                    paste(round(eigen(md$priorParameters$Lambda)$values, 6), collapse = ", ")))
      }
      
      cat(sprintf("     ✓ %s specification successful\n", spec_name))
    }
    
  }, error = function(e) {
    cat(sprintf("   ✗ Error in %s: %s\n", model_name, e$message))
  })
}

# Test all models
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
for (model in models) {
  debug_mvnormal_create(model)
}

# Test PosteriorParameters function
cat("\n--- Testing PosteriorParameters Function ---\n")
test_posterior_parameters <- function() {
  x <- matrix(rnorm(20), ncol = 2)
  
  tryCatch({
    # Create explicit parameters
    md <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = "EII"
    ))
    
    cat("Prior parameters:\n")
    cat(sprintf("   mu0: %s\n", paste(md$priorParameters$mu0, collapse = ", ")))
    cat(sprintf("   Lambda: %dx%d\n", nrow(md$priorParameters$Lambda), ncol(md$priorParameters$Lambda)))
    
    # Test PosteriorParameters
    post_params <- PosteriorParameters(md, x)
    
    cat("Posterior parameters:\n")
    cat(sprintf("   mu_n: %s\n", paste(round(post_params$mu_n, 3), collapse = ", ")))
    cat(sprintf("   t_n: %dx%d\n", nrow(post_params$t_n), ncol(post_params$t_n)))
    cat(sprintf("   nu_n: %g\n", post_params$nu_n))
    
    # Test rWishart with these parameters
    result <- rWishart(1, post_params$nu_n, post_params$t_n)
    cat("   ✓ rWishart successful\n")
    
  }, error = function(e) {
    cat(sprintf("   ✗ PosteriorParameters failed: %s\n", e$message))
  })
}

test_posterior_parameters()

# Test with different data
cat("\n--- Testing with Different Data Properties ---\n")
test_different_data <- function() {
  # Test with well-conditioned data
  data_cases <- list(
    normal = matrix(rnorm(20), ncol = 2),
    scaled = matrix(rnorm(20, sd = 0.1), ncol = 2),
    large = matrix(rnorm(100), ncol = 2),
    correlated = {
      x1 <- rnorm(10)
      x2 <- x1 + rnorm(10, sd = 0.1)
      cbind(x1, x2)
    }
  )
  
  for (case_name in names(data_cases)) {
    cat(sprintf("Testing %s data...\n", case_name))
    
    tryCatch({
      x <- data_cases[[case_name]]
      
      # Create with explicit parameters
      md <- MvnormalCreate(list(
        mu0 = c(0, 0),
        kappa0 = 1,
        nu = 3,
        Lambda = diag(2),
        covModel = "EII"
      ))
      
      # Test posterior parameters
      post_params <- PosteriorParameters(md, x)
      
      # Test rWishart
      result <- rWishart(1, post_params$nu_n, post_params$t_n)
      cat(sprintf("   ✓ %s data successful\n", case_name))
      
    }, error = function(e) {
      cat(sprintf("   ✗ %s data failed: %s\n", case_name, e$message))
    })
  }
}

test_different_data()

cat("\n=== PARAMETER CREATION DEBUG COMPLETE ===\n")