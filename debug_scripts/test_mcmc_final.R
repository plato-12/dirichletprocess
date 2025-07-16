devtools::load_all()
set_use_cpp(TRUE)

library(mvtnorm)
set.seed(42)
x <- matrix(rnorm(20), ncol = 2)

cat('Testing MCMC with updated Likelihood.mvnormal function...\n')
mdObj <- MvnormalCreate(list(
  mu0 = c(0, 0),
  kappa0 = 1,
  nu = 3,
  Lambda = diag(2),
  covModel = 'FULL'
))
dp <- DirichletProcessCreate(x, mdObj)
dp <- Initialise(dp)

cat('Before MCMC:\n')
ll_before <- LikelihoodDP(dp)
cat('LikelihoodDP success, length:', length(ll_before), '\n')

cat('Running 5 MCMC iterations...\n')
tryCatch({
  dp_after <- Fit(dp, 5)
  ll_after <- LikelihoodDP(dp_after)
  cat('SUCCESS: MCMC completed!\n')
  cat('Clusters after:', dp_after$numberClusters, '\n')
  cat('LikelihoodDP after MCMC, length:', length(ll_after), '\n')
}, error = function(e) {
  cat('ERROR:', e$message, '\n')
})