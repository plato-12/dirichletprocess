#!/usr/bin/env Rscript

library(dirichletprocess)
require(mvtnorm)

# Helper function for tolerance-based comparison
all_in_with_tolerance <- function(x, y, tolerance = 1e-10) {
  # For each element in x, check if there's at least one element in y that's close enough
  all(sapply(x, function(xi) any(abs(y - xi) < tolerance)))
}

# Test hierarchical mvnormal
dataTest <- list(rmvnorm(100, c(0,0), diag(2)), rmvnorm(100, c(1,1), diag(2)), rmvnorm(100, c(-1,-1), diag(2)), rmvnorm(100, c(2,2), diag(2)), rmvnorm(100, c(-2,-2), diag(2)))
dpobjlistTest <- DirichletProcessHierarchicalMvnormal2(dataTest)

cat("Before ClusterComponentUpdate:\n")
cat("Global parameters structure:\n")
cat("  Length:", length(dpobjlistTest$globalParameters), "\n")
if (length(dpobjlistTest$globalParameters) > 0) {
  cat("  Param 1 length:", length(dpobjlistTest$globalParameters[[1]]), "\n")
  cat("  Param 2 length:", length(dpobjlistTest$globalParameters[[2]]), "\n")
}

# Run update
dpobjlistTest <- ClusterComponentUpdate(dpobjlistTest)

cat("\nAfter ClusterComponentUpdate:\n")
cat("Global parameters structure:\n")
cat("  Length:", length(dpobjlistTest$globalParameters), "\n")
if (length(dpobjlistTest$globalParameters) > 0) {
  cat("  Param 1 length:", length(dpobjlistTest$globalParameters[[1]]), "\n")
  cat("  Param 2 length:", length(dpobjlistTest$globalParameters[[2]]), "\n")
}

# Check if global parameters contain the individual parameters
global_mu <- c(dpobjlistTest$globalParameters[[1]])
global_sigma <- c(dpobjlistTest$globalParameters[[2]])

cat("\nGlobal parameters:\n")
cat("  mu values (first 10):", paste(head(global_mu, 10), collapse=", "), "\n")
cat("  sigma values (first 10):", paste(head(global_sigma, 10), collapse=", "), "\n")

# Check all individual DPs for debugging
for (i in 1:length(dpobjlistTest$indDP)) {
  ind_mu <- c(dpobjlistTest$indDP[[i]]$clusterParameters[[1]])
  ind_sigma <- c(dpobjlistTest$indDP[[i]]$clusterParameters[[2]])
  
  cat("\nIndividual DP", i, ":\n")
  cat("  mu values:", paste(ind_mu, collapse=", "), "\n") 
  cat("  sigma values:", paste(ind_sigma, collapse=", "), "\n")
  
  # Test tolerance function
  mu_result <- all_in_with_tolerance(ind_mu, global_mu, tolerance = 2.0)
  sigma_result <- all_in_with_tolerance(ind_sigma, global_sigma, tolerance = 2.0)
  
  cat("  mu tolerance check:", mu_result, "\n")
  cat("  sigma tolerance check:", sigma_result, "\n")
  
  if (!sigma_result) {
    cat("  sigma differences (all values):\n")
    for (j in seq_along(ind_sigma)) {
      min_diff <- min(abs(global_sigma - ind_sigma[j]))
      cat("    ind_sigma[", j, "] =", ind_sigma[j], ", min global diff =", min_diff, "\n")
      if (min_diff >= 2.0) {
        cat("      *** EXCEEDS TOLERANCE ***\n")
      }
    }
  }
}