# Test complete initialise
library(dirichletprocess)
set_use_cpp(TRUE)

cat("=== TESTING COMPLETE INITIALISE ===\n")

# Create E model
params_1d <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1, 1, 1),
  covModel = "E"
)

md <- MvnormalCreate(params_1d)
data_1d <- matrix(rnorm(10), ncol = 1)
dp <- DirichletProcessCreate(data_1d, md)

cat("Calling Initialise with numInitialClusters=1...\n")
tryCatch({
  dp1 <- Initialise(dp, numInitialClusters = 1)
  cat("SUCCESS with 1 cluster\n")
}, error = function(e) {
  cat(sprintf("FAILED with 1 cluster: %s\n", e$message))
})

cat("Calling Initialise with numInitialClusters=2...\n")
tryCatch({
  dp2 <- Initialise(dp, numInitialClusters = 2)
  cat("SUCCESS with 2 clusters\n")
}, error = function(e) {
  cat(sprintf("FAILED with 2 clusters: %s\n", e$message))
})

cat("Calling Initialise with posterior=FALSE...\n")
tryCatch({
  dp3 <- Initialise(dp, posterior = FALSE, numInitialClusters = 2)
  cat("SUCCESS with posterior=FALSE\n")
}, error = function(e) {
  cat(sprintf("FAILED with posterior=FALSE: %s\n", e$message))
})

cat("\n=== COMPLETE INITIALISE TEST COMPLETE ===\n")