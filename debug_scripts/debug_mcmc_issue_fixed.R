devtools::load_all()
set_use_cpp(FALSE)  # Start with R implementation

library(mvtnorm)
set.seed(42)
x <- matrix(rnorm(20), ncol = 2)

cat('=== Testing with R implementation first ===\n')
mdObj <- MvnormalCreate(list(
  mu0 = c(0, 0),
  kappa0 = 1,
  nu = 3,
  Lambda = diag(2),
  covModel = "FULL"
))
dp <- DirichletProcessCreate(x, mdObj)
dp <- Initialise(dp)

cat('R implementation - Before MCMC:\n')
cat('Number of clusters:', dp$numberClusters, '\n')
ll_before <- LikelihoodDP(dp)
cat('LikelihoodDP success, length:', length(ll_before), '\n')

cat('R implementation - Running 1 MCMC iteration:\n')
tryCatch({
  dp_r <- Fit(dp, 1)
  ll_after <- LikelihoodDP(dp_r)
  cat('MCMC successful, clusters after:', dp_r$numberClusters, '\n')
}, error = function(e) {
  cat('MCMC failed with R:', e$message, '\n')
})

cat('\n=== Now testing with C++ implementation ===\n')
set_use_cpp(TRUE)

mdObj_cpp <- MvnormalCreate(list(
  mu0 = c(0, 0),
  kappa0 = 1,
  nu = 3,
  Lambda = diag(2),
  covModel = "FULL"
))
dp_cpp <- DirichletProcessCreate(x, mdObj_cpp)
dp_cpp <- Initialise(dp_cpp)

cat('C++ implementation - Before MCMC:\n')
cat('Number of clusters:', dp_cpp$numberClusters, '\n')
ll_before_cpp <- LikelihoodDP(dp_cpp)
cat('LikelihoodDP success, length:', length(ll_before_cpp), '\n')

cat('C++ implementation - Running 1 MCMC iteration:\n')
tryCatch({
  dp_cpp_after <- Fit(dp_cpp, 1)
  ll_after_cpp <- LikelihoodDP(dp_cpp_after)
  cat('MCMC successful, clusters after:', dp_cpp_after$numberClusters, '\n')
}, error = function(e) {
  cat('MCMC failed with C++:', e$message, '\n')
})