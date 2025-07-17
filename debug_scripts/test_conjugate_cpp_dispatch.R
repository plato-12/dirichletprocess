#!/usr/bin/env Rscript

# Test conjugate C++ dispatch for normal distribution
library(dirichletprocess)

cat("Testing C++ dispatch for conjugate normal distribution\n")
cat("===================================================\n")

# Enable C++ 
set_use_cpp(TRUE)
cat("C++ enabled:", using_cpp(), "\n")

# Create normal distribution DP
y <- rnorm(50)
dp <- DirichletProcessGaussian(y)

cat("DP classes:", paste(class(dp), collapse=" "), "\n")
cat("Mixing distribution classes:", paste(class(dp$mixingDistribution), collapse=" "), "\n")

# Check S3 method dispatch availability
cat("\nChecking S3 method dispatch:\n")
cat("ClusterComponentUpdate methods:\n")
methods_cc <- methods("ClusterComponentUpdate")
print(methods_cc)

cat("\nClusterParameterUpdate methods:\n")
methods_cp <- methods("ClusterParameterUpdate")
print(methods_cp)

# Test if C++ methods are actually being called
cat("\nTesting method dispatch:\n")

# Temporarily add a class to force C++ dispatch
cat("Adding 'conjugate.cpp' class to test C++ dispatch...\n")
test_dp <- dp
class(test_dp) <- c(class(test_dp), "conjugate.cpp")
cat("Modified classes:", paste(class(test_dp), collapse=" "), "\n")

# Test C++ method dispatch
tryCatch({
  cat("Testing ClusterComponentUpdate with cpp class...\n")
  test_dp <- ClusterComponentUpdate(test_dp)
  cat("SUCCESS: ClusterComponentUpdate with cpp class worked\n")
}, error = function(e) {
  cat("ERROR in ClusterComponentUpdate with cpp class:", e$message, "\n")
})

tryCatch({
  cat("Testing ClusterParameterUpdate with cpp class...\n")
  test_dp <- ClusterParameterUpdate(test_dp)
  cat("SUCCESS: ClusterParameterUpdate with cpp class worked\n")
}, error = function(e) {
  cat("ERROR in ClusterParameterUpdate with cpp class:", e$message, "\n")
})

# Test manual MCMC with potential C++ dispatch
cat("\nTesting manual MCMC loop with original classes:\n")
samples <- list()
for(s in seq_len(5)){
  dp <- ClusterComponentUpdate(dp)
  dp <- ClusterParameterUpdate(dp)
  if(s %% 3 == 0) {
    dp <- UpdateAlpha(dp)
  }
  samples[[s]] <- list()
  samples[[s]]$phi <- dp$clusterParameters
  samples[[s]]$weights <- dp$pointsPerCluster / dp$n
  
  cat("Iteration", s, "- clusters:", dp$numberClusters, "\n")
}

cat("\nManual MCMC completed successfully!\n")
cat("Final clusters:", dp$numberClusters, "\n")

# Check which C++ functions are actually available
cat("\nDirect C++ function availability:\n")
cat("conjugate_cluster_component_update_cpp in namespace:", 
    exists("conjugate_cluster_component_update_cpp", where=getNamespace("dirichletprocess")), "\n")
cat("conjugate_cluster_parameter_update_cpp in namespace:",
    exists("conjugate_cluster_parameter_update_cpp", where=getNamespace("dirichletprocess")), "\n")

cat("\nTest completed!\n")