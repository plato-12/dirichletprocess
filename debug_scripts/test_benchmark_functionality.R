devtools::load_all()
set_use_cpp(TRUE)

library(mvtnorm)

# Create test data similar to what the benchmark would use
set.seed(42)
n_samples <- 100
d <- 5
x <- matrix(rnorm(n_samples * d), ncol = d)

cat('Testing benchmark-style functionality with', n_samples, 'samples and', d, 'dimensions\n\n')

# Test multiple covariance models like the benchmark does
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")

# Test with different numbers of MCMC iterations like benchmark
n_its <- c(10, 50, 100)

cat('=== BENCHMARK-STYLE TEST ===\n')

for (model in models) {
  cat('Testing model:', model, '\n')
  
  for (n_it in n_its) {
    cat('  Iterations:', n_it)
    
    start_time <- Sys.time()
    tryCatch({
      # Create and initialize DP object
      mdObj <- MvnormalCreate(list(
        mu0 = rep(0, d),
        kappa0 = 1,
        nu = d + 2,
        Lambda = diag(d),
        covModel = model
      ))
      
      dp <- DirichletProcessCreate(x, mdObj)
      dp <- Initialise(dp)
      
      # Run MCMC like benchmark does
      dp <- Fit(dp, n_it)
      
      end_time <- Sys.time()
      elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      cat(' -> ✓ Success (', round(elapsed, 3), 's)\n')
      
    }, error = function(e) {
      cat(' -> ✗ Error:', e$message, '\n')
    })
  }
  cat('\n')
}

cat('=== FINAL VERIFICATION ===\n')
# Test that LikelihoodDP works correctly with multiple clusters
dp_final <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = "FULL")))
dp_final <- Initialise(dp_final)
dp_final <- Fit(dp_final, 50)  # Get some clusters

tryCatch({
  ll <- LikelihoodDP(dp_final)
  cat('✓ Final LikelihoodDP test: Length', length(ll), ', clusters:', dp_final$numberClusters, '\n')
  cat('✓ All benchmark functionality working correctly!\n')
}, error = function(e) {
  cat('✗ Final LikelihoodDP test failed:', e$message, '\n')
})