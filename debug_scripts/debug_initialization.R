devtools::load_all()
set_use_cpp(TRUE)

library(mvtnorm)
set.seed(42)
x <- matrix(rnorm(20), ncol=2)

# Test step by step
cat('Step 1: Creating mvnormal mixing distribution...\n')
mdObj <- MvnormalCreate(list(
  mu0 = c(0, 0),
  kappa0 = 1,
  nu = 3,
  Lambda = diag(2),
  covModel = 'FULL'
))

cat('Step 2: Creating DirichletProcess object...\n')
dp <- DirichletProcessCreate(x, mdObj)

cat('Step 3: Initializing DirichletProcess object...\n')
dp <- Initialise(dp)

cat('Step 4: Check dp structure after initialization...\n')
cat('Number of clusters:', dp$numberClusters, '\n')
cat('Data dimensions:', dim(dp$data), '\n')

cat('Step 5: Check cluster parameters after initialization...\n')
cat('Cluster parameters names:', names(dp$clusterParameters), '\n')
if ('mu' %in% names(dp$clusterParameters)) {
  cat('mu dimensions:', dim(dp$clusterParameters$mu), '\n')
  cat('mu values:', dp$clusterParameters$mu[,,1], '\n')
}
if ('sig' %in% names(dp$clusterParameters)) {
  cat('sig dimensions:', dim(dp$clusterParameters$sig), '\n')
}

cat('Step 6: Test direct likelihood call after initialization...\n')
tryCatch({
  lik_result <- Likelihood(mdObj, x[1, , drop=FALSE], dp$clusterParameters)
  cat('Direct Likelihood result:', lik_result, '\n')
}, error = function(e) {
  cat('Direct Likelihood error:', e$message, '\n')
})

cat('Step 7: Test LikelihoodDP after initialization...\n')
tryCatch({
  ll <- LikelihoodDP(dp)
  cat('LikelihoodDP success, length:', length(ll), '\n')
  cat('First few values:', ll[1:min(5, length(ll))], '\n')
}, error = function(e) {
  cat('LikelihoodDP error:', e$message, '\n')
})