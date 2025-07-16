# Test script to validate constrained covariance fixes
# Tests both cluster expansion and dimension access issues

library(dirichletprocess)
set.seed(123)

# Test cluster expansion with large datasets  
test_cluster_expansion <- function() {
  cat("Testing cluster expansion for constrained models...\n")
  
  constrained_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  
  results <- list()
  
  for (model in constrained_models) {
    cat(sprintf("Testing %s model...\n", model))
    
    # Create larger dataset to force cluster expansion
    n <- 100
    x <- matrix(rnorm(n * 2), ncol = 2)
    
    tryCatch({
      # Create DP object
      dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model)))
      dp <- Initialise(dp)
      
      # Run MCMC for multiple iterations to test expansion
      dp <- Fit(dp, its = 5, progressBar = FALSE)
      
      results[[model]] <- list(
        success = TRUE,
        final_clusters = dp$numberClusters,
        message = sprintf("Success: %d clusters created", dp$numberClusters)
      )
      cat(sprintf("  ✓ %s: %d clusters\n", model, dp$numberClusters))
      
    }, error = function(e) {
      results[[model]] <- list(
        success = FALSE,
        error = e$message,
        message = sprintf("Failed: %s", e$message)
      )
      cat(sprintf("  ✗ %s: %s\n", model, e$message))
    })
  }
  
  return(results)
}

# Test complete MCMC pipeline with various scenarios
test_mcmc_pipeline <- function() {
  cat("\nTesting complete MCMC pipeline...\n")
  
  constrained_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
  
  results <- list()
  
  for (model in constrained_models) {
    cat(sprintf("Testing %s model with 10 iterations...\n", model))
    
    # Create test data
    x <- matrix(rnorm(50), ncol = 2)
    
    tryCatch({
      # Create DP object
      dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model)))
      dp <- Initialise(dp)
      
      # Run extended MCMC
      dp <- Fit(dp, its = 10, progressBar = FALSE)
      
      results[[model]] <- list(
        success = TRUE,
        final_clusters = dp$numberClusters,
        iterations_completed = 10,
        message = sprintf("Success: %d clusters after 10 iterations", dp$numberClusters)
      )
      cat(sprintf("  ✓ %s: %d clusters after 10 iterations\n", model, dp$numberClusters))
      
    }, error = function(e) {
      results[[model]] <- list(
        success = FALSE,
        error = e$message,
        message = sprintf("Failed: %s", e$message)
      )
      cat(sprintf("  ✗ %s: %s\n", model, e$message))
    })
  }
  
  return(results)
}

# Test individual components
test_individual_components <- function() {
  cat("\nTesting individual MCMC components...\n")
  
  model <- "EII"  # Test with one constrained model
  x <- matrix(rnorm(20), ncol = 2)
  
  tryCatch({
    dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = model)))
    dp <- Initialise(dp)
    
    cat("  Testing ClusterComponentUpdate...\n")
    dp <- ClusterComponentUpdate(dp)
    cat("  ✓ ClusterComponentUpdate passed\n")
    
    cat("  Testing ClusterParameterUpdate...\n")
    dp <- ClusterParameterUpdate(dp)
    cat("  ✓ ClusterParameterUpdate passed\n")
    
    cat("  Testing UpdateAlpha...\n")
    dp <- UpdateAlpha(dp)
    cat("  ✓ UpdateAlpha passed\n")
    
    return(list(success = TRUE, message = "All components working"))
    
  }, error = function(e) {
    cat(sprintf("  ✗ Component test failed: %s\n", e$message))
    return(list(success = FALSE, error = e$message))
  })
}

# Main test execution
cat("=== Testing Constrained Covariance Fixes ===\n")

# Test 1: Cluster expansion
expansion_results <- test_cluster_expansion()

# Test 2: Complete MCMC pipeline
pipeline_results <- test_mcmc_pipeline()

# Test 3: Individual components
component_results <- test_individual_components()

# Summary
cat("\n=== SUMMARY ===\n")
cat(sprintf("Cluster expansion tests: %d/%d passed\n", 
            sum(sapply(expansion_results, function(x) x$success)), 
            length(expansion_results)))

cat(sprintf("MCMC pipeline tests: %d/%d passed\n", 
            sum(sapply(pipeline_results, function(x) x$success)), 
            length(pipeline_results)))

cat(sprintf("Individual component tests: %s\n", 
            if(component_results$success) "PASSED" else "FAILED"))

# Report any failures
failures <- c()
for (model in names(expansion_results)) {
  if (!expansion_results[[model]]$success) {
    failures <- c(failures, sprintf("Expansion %s: %s", model, expansion_results[[model]]$error))
  }
}
for (model in names(pipeline_results)) {
  if (!pipeline_results[[model]]$success) {
    failures <- c(failures, sprintf("Pipeline %s: %s", model, pipeline_results[[model]]$error))
  }
}
if (!component_results$success) {
  failures <- c(failures, sprintf("Components: %s", component_results$error))
}

if (length(failures) > 0) {
  cat("\nFAILURES:\n")
  for (failure in failures) {
    cat(sprintf("  - %s\n", failure))
  }
} else {
  cat("\n✓ ALL TESTS PASSED!\n")
}