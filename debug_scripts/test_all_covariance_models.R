devtools::load_all()
set_use_cpp(TRUE)

library(mvtnorm)
set.seed(42)
x <- matrix(rnorm(20), ncol=2)

# Test all covariance models
cat('Testing all covariance models with C++ enabled...\n\n')
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")

results <- list()

for (model in models) {
  cat('Testing model:', model, '\n')
  tryCatch({
    # Create mixing distribution
    mdObj <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = model
    ))
    
    # Create and initialize DP object
    dp <- DirichletProcessCreate(x, mdObj)
    dp <- Initialise(dp)
    
    # Test likelihood
    ll <- LikelihoodDP(dp)
    
    cat('  ✓ Success: Length', length(ll), ', sample values:', ll[1:min(3, length(ll))], '\n')
    results[[model]] <- list(success = TRUE, length = length(ll), values = ll[1:3])
    
  }, error = function(e) {
    cat('  ✗ Error:', e$message, '\n')
    results[[model]] <- list(success = FALSE, error = e$message)
  })
  cat('\n')
}

# Summary
cat('\n=== SUMMARY ===\n')
successes <- sum(sapply(results, function(x) x$success))
cat('Successful models:', successes, 'out of', length(models), '\n')

if (successes > 0) {
  cat('Working models:', paste(names(results)[sapply(results, function(x) x$success)], collapse = ', '), '\n')
}

if (successes < length(models)) {
  failed <- names(results)[!sapply(results, function(x) x$success)]
  cat('Failed models:', paste(failed, collapse = ', '), '\n')
}

cat('\nC++ mvnormal implementation status: ')
if (successes == length(models)) {
  cat('✓ ALL WORKING\n')
} else if (successes > 0) {
  cat('⚠ PARTIALLY WORKING\n')
} else {
  cat('✗ NOT WORKING\n')
}