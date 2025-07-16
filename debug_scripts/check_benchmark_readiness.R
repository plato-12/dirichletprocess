# Check if benchmark is ready to run
library(dirichletprocess)
set.seed(123)

cat("=== CHECKING BENCHMARK READINESS ===\n")

# Test 1: Can we run the atime benchmark function?
test_atime_benchmark <- function() {
  cat("Testing atime benchmark execution...\n")
  
  tryCatch({
    # Source the benchmark file
    source("benchmark/atime/benchmark-covariance-models-comprehensive.R")
    
    # Try to run the atime benchmark
    cat("Attempting to run atime benchmark...\n")
    
    # This should work if everything is fixed
    atime_results <- run_atime_benchmark()
    
    cat("✓ Atime benchmark completed successfully!\n")
    cat(sprintf("Results: %d benchmark comparisons\n", length(atime_results)))
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("✗ Atime benchmark failed: %s\n", e$message))
    cat("Full error details:\n")
    print(e)
    return(FALSE)
  })
}

# Test 2: Are individual models working?
test_individual_models <- function() {
  cat("\nTesting individual model functionality...\n")
  
  models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
  x <- matrix(rnorm(40), ncol = 2)
  
  success_count <- 0
  
  for (model in models) {
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
      dp <- Fit(dp, 10, progressBar = FALSE)
      
      success_count <- success_count + 1
      cat(sprintf("  ✓ %s: SUCCESS\n", model))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s: FAILED - %s\n", model, e$message))
    })
  }
  
  cat(sprintf("Individual models: %d/%d working\n", success_count, length(models)))
  return(success_count == length(models))
}

# Test 3: Check data generation
test_data_generation <- function() {
  cat("\nTesting data generation...\n")
  
  tryCatch({
    # Test the data generation function
    source("benchmark/atime/benchmark-covariance-models-comprehensive.R")
    
    # Test different dimensions
    test_data_1d <- generate_benchmark_data(50, 1)
    test_data_2d <- generate_benchmark_data(50, 2)
    test_data_5d <- generate_benchmark_data(50, 5)
    
    cat(sprintf("  ✓ 1D data: %s\n", paste(dim(test_data_1d), collapse = "x")))
    cat(sprintf("  ✓ 2D data: %s\n", paste(dim(test_data_2d), collapse = "x")))
    cat(sprintf("  ✓ 5D data: %s\n", paste(dim(test_data_5d), collapse = "x")))
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("  ✗ Data generation failed: %s\n", e$message))
    return(FALSE)
  })
}

# Run all tests
individual_ok <- test_individual_models()
data_gen_ok <- test_data_generation()
benchmark_ok <- test_atime_benchmark()

cat("\n=== READINESS ASSESSMENT ===\n")
cat(sprintf("Individual models: %s\n", if (individual_ok) "✓ READY" else "✗ NOT READY"))
cat(sprintf("Data generation: %s\n", if (data_gen_ok) "✓ READY" else "✗ NOT READY"))
cat(sprintf("Atime benchmark: %s\n", if (benchmark_ok) "✓ READY" else "✗ NOT READY"))

if (individual_ok && data_gen_ok && benchmark_ok) {
  cat("\n🎉 ALL SYSTEMS READY! You can run atime_results <- run_atime_benchmark()\n")
} else {
  cat("\n⚠️  ISSUES REMAIN - Do not run atime benchmark yet\n")
}

cat("\n=== READINESS CHECK COMPLETE ===\n")