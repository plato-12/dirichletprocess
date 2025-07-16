# Final Verification Test
# ======================
# 
# This test verifies that the MVNormal C++ parameter mismatch issue is fully resolved.
# Tests the exact scenario mentioned in the analysis document:
# - Multiple covariance models
# - C++ enabled
# - MCMC with multiple iterations
# - LikelihoodDP functionality during MCMC

devtools::load_all()
set_use_cpp(TRUE)

library(mvtnorm)

cat('=== FINAL VERIFICATION: MVNormal C++ Parameter Mismatch Fix ===\n\n')

# Create test data similar to benchmark conditions
set.seed(42)
n_samples <- 100
d <- 5
x <- matrix(rnorm(n_samples * d), ncol = d)

cat('Test configuration:\n')
cat('- Samples:', n_samples, '\n')
cat('- Dimensions:', d, '\n')
cat('- C++ enabled:', using_cpp(), '\n')
cat('- C++ samplers enabled:', using_cpp_samplers(), '\n\n')

# Test all 7 covariance models with MCMC
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
results <- list()

for (model in models) {
  cat('Testing', model, 'covariance model:\n')
  
  start_time <- Sys.time()
  success <- tryCatch({
    # Create DP object
    mdObj <- MvnormalCreate(list(
      mu0 = rep(0, d),
      kappa0 = 1,
      nu = d + 2,
      Lambda = diag(d),
      covModel = model
    ))
    
    dp <- DirichletProcessCreate(x, mdObj)
    dp <- Initialise(dp)
    
    # Test initial LikelihoodDP
    ll_before <- LikelihoodDP(dp)
    cat('  ✓ Initial LikelihoodDP: length', length(ll_before), '\n')
    
    # Run MCMC iterations
    dp_fitted <- Fit(dp, 20)  # 20 iterations to ensure cluster changes
    
    # Test final LikelihoodDP  
    ll_after <- LikelihoodDP(dp_fitted)
    cat('  ✓ MCMC completed: clusters', dp_fitted$numberClusters, '\n')
    cat('  ✓ Final LikelihoodDP: length', length(ll_after), '\n')
    
    # Verify results make sense
    if (length(ll_before) != length(ll_after) || length(ll_after) != n_samples) {
      stop("LikelihoodDP length mismatch")
    }
    
    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cat('  ✓ Time:', round(elapsed, 2), 'seconds\n')
    
    TRUE
  }, error = function(e) {
    cat('  ✗ ERROR:', e$message, '\n')
    FALSE
  })
  
  results[[model]] <- success
  cat('\n')
}

# Summary
cat('=== SUMMARY ===\n')
successes <- sum(unlist(results))
total <- length(results)

cat('Successful models:', successes, 'out of', total, '\n')

if (successes > 0) {
  working <- names(results)[unlist(results)]
  cat('Working models:', paste(working, collapse = ', '), '\n')
}

if (successes < total) {
  failed <- names(results)[!unlist(results)]
  cat('Failed models:', paste(failed, collapse = ', '), '\n')
}

# Final verdict
cat('\n=== FINAL VERDICT ===\n')
if (successes == total) {
  cat('🎉 SUCCESS: All MVNormal C++ parameter mismatch issues RESOLVED!\n')
  cat('✅ All', total, 'covariance models work correctly with C++ enabled\n')
  cat('✅ MCMC functionality fully working\n') 
  cat('✅ LikelihoodDP functionality fully working\n')
  cat('✅ Benchmark functionality should now work correctly\n')
} else if (successes > 0) {
  cat('⚠️  PARTIAL SUCCESS:', successes, 'out of', total, 'models working\n')
  cat('❌ Some issues remain for:', paste(names(results)[!unlist(results)], collapse = ', '), '\n')
} else {
  cat('❌ FAILURE: MVNormal C++ parameter mismatch issues NOT resolved\n')
}

cat('\nOriginal issue: "values must be length 2, but FUN(X[[1]]) result is length 1"\n')
cat('Status: Should be RESOLVED if all tests passed\n')