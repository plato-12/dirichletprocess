# Debug the exact benchmark issue
library(dirichletprocess)

# Disable C++ to test R implementation
set_use_cpp(FALSE)

# Create test data similar to benchmark
set.seed(42)
test_data <- matrix(rnorm(50 * 5), ncol = 5)  # 50 samples, 5 features

# Create prior parameters for FULL model
prior_params <- list(
  mu0 = rep(0, 5),
  kappa0 = 1,
  nu = 7,
  Lambda = diag(5),
  covModel = "FULL"
)

cat("=== BENCHMARK REPLICATION ===\n")
tryCatch({
  # Create mixing distribution and Dirichlet process
  md <- MvnormalCreate(prior_params)
  dp <- DirichletProcessCreate(test_data, md)
  dp <- Initialise(dp, numInitialClusters = 2)
  
  cat("After initialization:\n")
  cat("  numberClusters:", dp$numberClusters, "\n")
  cat("  mu dimensions:", dim(dp$clusterParameters$mu), "\n")
  cat("  sig dimensions:", dim(dp$clusterParameters$sig), "\n")
  
  # Test LikelihoodDP after initialization
  cat("Testing LikelihoodDP after initialization...\n")
  lik_dp <- LikelihoodDP(dp)
  cat("  LikelihoodDP result dimensions:", dim(lik_dp), "\n")
  
  # Run just one iteration of Fit to see where it fails
  cat("Running one iteration of Fit...\n")
  dp <- Fit(dp, 1, progressBar = FALSE)
  
  cat("After 1 iteration:\n")
  cat("  numberClusters:", dp$numberClusters, "\n")
  cat("  mu dimensions:", dim(dp$clusterParameters$mu), "\n")
  cat("  sig dimensions:", dim(dp$clusterParameters$sig), "\n")
  
  # Test LikelihoodDP after one iteration
  cat("Testing LikelihoodDP after one iteration...\n")
  lik_dp2 <- LikelihoodDP(dp)
  cat("  LikelihoodDP result dimensions:", dim(lik_dp2), "\n")
  
  cat("SUCCESS: Benchmark replication worked!\n")
  
}, error = function(e) {
  cat("ERROR in benchmark replication:", e$message, "\n")
  cat("Call stack:\n")
  traceback()
})

cat("\n=== TESTING DIFFERENT COVARIANCE MODELS ===\n")

# Test different covariance models
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")

for (model in models) {
  cat(sprintf("Testing model: %s\n", model))
  
  # Create prior parameters for this model
  prior_params <- list(
    mu0 = rep(0, 5),
    kappa0 = 1,
    nu = 7,
    Lambda = diag(5),
    covModel = model
  )
  
  tryCatch({
    # Create mixing distribution and Dirichlet process
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(test_data, md)
    dp <- Initialise(dp, numInitialClusters = 2)
    
    # Test LikelihoodDP
    lik_dp <- LikelihoodDP(dp)
    
    # Run one iteration of Fit
    dp <- Fit(dp, 1, progressBar = FALSE)
    
    cat(sprintf("  %s: SUCCESS\n", model))
    
  }, error = function(e) {
    cat(sprintf("  %s: ERROR - %s\n", model, e$message))
  })
}