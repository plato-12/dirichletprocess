# Test bypass predictive
library(dirichletprocess)
set_use_cpp(FALSE)  # Turn off C++

cat("=== TESTING BYPASS PREDICTIVE ===\n")

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

# Manual Initialise.conjugate without InitialisePredictive
cat("Manual Initialise.conjugate without InitialisePredictive...\n")
tryCatch({
  numInitialClusters <- 2
  dp$clusterLabels <- rep_len(seq_len(numInitialClusters), length.out = dp$n)
  dp$numberClusters <- numInitialClusters
  dp$pointsPerCluster <- vapply(seq_len(numInitialClusters), function(x) sum(dp$clusterLabels == x), numeric(1))
  
  # This should work
  dp$clusterParameters <- PriorDraw(dp$mixingDistribution, numInitialClusters)
  
  cat("Manual initialise without predictive SUCCESS\n")
  cat(sprintf("Clusters: %d\n", dp$numberClusters))
  
}, error = function(e) {
  cat(sprintf("Manual initialise FAILED: %s\n", e$message))
})

# Now test just the InitialisePredictive call
cat("\nTesting InitialisePredictive call...\n")
tryCatch({
  # This should also work since Predictive worked
  dp$predictiveArray <- Predictive(dp$mixingDistribution, dp$data)
  cat("InitialisePredictive SUCCESS\n")
  
}, error = function(e) {
  cat(sprintf("InitialisePredictive FAILED: %s\n", e$message))
})

cat("\n=== BYPASS PREDICTIVE TEST COMPLETE ===\n")