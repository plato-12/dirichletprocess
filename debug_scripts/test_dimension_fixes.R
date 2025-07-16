# Test dimension access fixes specifically
# Focus on verifying the fixes work without running full MCMC

library(dirichletprocess)
set.seed(123)

# Test basic initialization and parameter access
test_basic_initialization <- function() {
  cat("Testing basic initialization and parameter access...\n")
  
  constrained_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  
  for (model in constrained_models) {
    cat(sprintf("Testing %s model...\n", model))
    
    tryCatch({
      # Create simple test data
      x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
      
      # Create DP object
      dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model)))
      dp <- Initialise(dp)
      
      # Test parameter dimensions
      param1_dims <- dim(dp$clusterParameters[[1]])
      param2_dims <- dim(dp$clusterParameters[[2]])
      
      cat(sprintf("  Parameter 1 dimensions: %s\n", paste(param1_dims, collapse = "x")))
      cat(sprintf("  Parameter 2 dimensions: %s\n", paste(param2_dims, collapse = "x")))
      
      # Test parameter access patterns
      for (i in 1:dp$numberClusters) {
        # Test dimension-aware parameter access
        if (length(param1_dims) == 3) {
          mu_i <- dp$clusterParameters[[1]][, , i]
        } else if (length(param1_dims) == 2) {
          mu_i <- dp$clusterParameters[[1]][, i]
        } else {
          mu_i <- dp$clusterParameters[[1]][i]
        }
        
        if (length(param2_dims) == 3) {
          sig_i <- dp$clusterParameters[[2]][, , i]
        } else if (length(param2_dims) == 2) {
          sig_i <- dp$clusterParameters[[2]][, i]
        } else {
          sig_i <- dp$clusterParameters[[2]][i]
        }
        
        cat(sprintf("  Cluster %d: mu shape = %s, sig shape = %s\n", 
                    i, paste(dim(mu_i), collapse = "x"), paste(dim(sig_i), collapse = "x")))
      }
      
      cat(sprintf("  ✓ %s: Initialization successful\n", model))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s: %s\n", model, e$message))
    })
  }
}

# Test likelihood calculations with constrained models
test_likelihood_calculations <- function() {
  cat("\nTesting likelihood calculations...\n")
  
  constrained_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  
  for (model in constrained_models) {
    cat(sprintf("Testing %s model likelihood...\n", model))
    
    tryCatch({
      # Create simple test data
      x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
      
      # Create DP object
      dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model)))
      dp <- Initialise(dp)
      
      # Test likelihood calculation
      theta <- dp$clusterParameters
      x_test <- matrix(c(1.5, 2.5), nrow = 1)
      
      # This should work with the dimension-aware fixes
      likelihood_result <- Likelihood(dp$mixingDistribution, x_test, theta)
      
      cat(sprintf("  ✓ %s: Likelihood calculation successful\n", model))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s: %s\n", model, e$message))
    })
  }
}

# Test parameter updates
test_parameter_updates <- function() {
  cat("\nTesting parameter updates...\n")
  
  # Test with one model to avoid LAPACK issues
  model <- "EII"
  
  tryCatch({
    # Create simple test data
    x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
    
    # Create DP object
    dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model)))
    dp <- Initialise(dp)
    
    # Test PosteriorDraw
    x_test <- matrix(c(1.5, 2.5), nrow = 1)
    post_draw <- PosteriorDraw(dp$mixingDistribution, x_test)
    
    cat(sprintf("  Posterior draw successful: %d parameters\n", length(post_draw)))
    
    # Test manual parameter assignment with dimension awareness
    param1_dims <- dim(dp$clusterParameters[[1]])
    param2_dims <- dim(dp$clusterParameters[[2]])
    
    if (length(param1_dims) == 3) {
      dp$clusterParameters[[1]][, , 1] <- post_draw[[1]]
    } else if (length(param1_dims) == 2) {
      dp$clusterParameters[[1]][, 1] <- post_draw[[1]]
    } else {
      dp$clusterParameters[[1]][1] <- post_draw[[1]]
    }
    
    if (length(param2_dims) == 3) {
      dp$clusterParameters[[2]][, , 1] <- post_draw[[2]]
    } else if (length(param2_dims) == 2) {
      dp$clusterParameters[[2]][, 1] <- post_draw[[2]]
    } else {
      dp$clusterParameters[[2]][1] <- post_draw[[2]]
    }
    
    cat(sprintf("  ✓ Parameter assignment successful\n"))
    
  }, error = function(e) {
    cat(sprintf("  ✗ Parameter update failed: %s\n", e$message))
  })
}

# Main test execution
cat("=== Testing Dimension Access Fixes ===\n")

# Test 1: Basic initialization
test_basic_initialization()

# Test 2: Likelihood calculations
test_likelihood_calculations()

# Test 3: Parameter updates
test_parameter_updates()

cat("\n=== Test Complete ===\n")