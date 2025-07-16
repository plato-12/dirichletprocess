devtools::load_all()
set_use_cpp(TRUE)

# Test for FULL covariance model
library(mvtnorm)
set.seed(42)
x <- matrix(rnorm(20), ncol=2)

# Create a DP object using MvnormalCreate
cat('Creating mvnormal mixing distribution with FULL covariance...\n')
mdObj <- MvnormalCreate(list(
  mu0 = c(0, 0),
  kappa0 = 1,
  nu = 3,
  Lambda = diag(2),
  covModel = 'FULL'
))

dp <- DirichletProcessCreate(x, mdObj)

# Test basic likelihood functionality
cat('Testing basic mvnormal likelihood with C++...\n')
tryCatch({
  ll <- LikelihoodDP(dp)
  cat('Success: LikelihoodDP returned vector of length', length(ll), '\n')
  cat('Likelihood values:', ll[1:min(5, length(ll))], '\n')
}, error = function(e) {
  cat('Error:', e$message, '\n')
})

# Test different covariance models
cat('\nTesting different covariance models...\n')
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")

for (model in models) {
  cat('Testing model:', model, '\n')
  tryCatch({
    mdObj <- MvnormalCreate(list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2),
      covModel = model
    ))
    
    dp <- DirichletProcessCreate(x, mdObj)
    ll <- LikelihoodDP(dp)
    cat('  Success: Length', length(ll), 'values:', ll[1:min(3, length(ll))], '\n')
  }, error = function(e) {
    cat('  Error:', e$message, '\n')
  })
}