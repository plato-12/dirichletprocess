#!/usr/bin/env Rscript

# Test manual MCMC loop with DirichletProcessGaussian (normal distribution)
library(dirichletprocess)

cat("Testing manual MCMC loop with DirichletProcessGaussian\n")
cat("=====================================================\n")

# Enable C++ if available
set_use_cpp(TRUE)
cat("C++ enabled:", using_cpp(), "\n")

# Create sample data and Dirichlet process
y <- rnorm(50)
dp <- DirichletProcessGaussian(y)

cat("Initial dp classes:", paste(class(dp), collapse=" "), "\n")
cat("Mixing distribution classes:", paste(class(dp$mixingDistribution), collapse=" "), "\n") 
cat("Initial clusters:", dp$numberClusters, "\n")

# Test individual MCMC functions
cat("\n1. Testing ClusterComponentUpdate...\n")
dp <- ClusterComponentUpdate(dp)
cat("   After ClusterComponentUpdate - clusters:", dp$numberClusters, "\n")

cat("2. Testing ClusterParameterUpdate...\n")
dp <- ClusterParameterUpdate(dp)  
cat("   After ClusterParameterUpdate - clusters:", dp$numberClusters, "\n")

cat("3. Testing UpdateAlpha...\n")
dp <- UpdateAlpha(dp)
cat("   After UpdateAlpha - alpha:", dp$alpha, "\n")

cat("4. Testing weights calculation...\n")
weights <- dp$pointsPerCluster / dp$n
cat("   Weights calculated successfully, length:", length(weights), "\n")
cat("   Weights sum:", sum(weights), "\n")

# Test the exact manual MCMC loop from vignette
cat("\n5. Testing complete manual MCMC loop (vignette style)...\n")
samples <- list()
for(s in seq_len(10)){
  dp <- ClusterComponentUpdate(dp)
  dp <- ClusterParameterUpdate(dp)
  if(s %% 5 == 0) {
    dp <- UpdateAlpha(dp)
  }
  samples[[s]] <- list()
  samples[[s]]$phi <- dp$clusterParameters
  samples[[s]]$weights <- dp$pointsPerCluster / dp$n
  
  cat("   Iteration", s, "- clusters:", dp$numberClusters, 
      "alpha:", round(dp$alpha, 3), "\n")
}

cat("\nManual MCMC loop completed successfully!\n")
cat("Number of samples:", length(samples), "\n")
cat("Final number of clusters:", dp$numberClusters, "\n")
cat("Final alpha:", dp$alpha, "\n")

# Check C++ function availability for normal distribution
cat("\nC++ Function Availability Check:\n")
cat("conjugate_cluster_component_update_cpp:", exists("conjugate_cluster_component_update_cpp"), "\n")
cat("conjugate_cluster_parameter_update_cpp:", exists("conjugate_cluster_parameter_update_cpp"), "\n")
cat("update_concentration function:", exists("update_concentration"), "\n")

cat("\nTest completed successfully!\n")