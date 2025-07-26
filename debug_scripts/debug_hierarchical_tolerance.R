#!/usr/bin/env Rscript

library(dirichletprocess)

# Helper function for tolerance-based comparison
all_in_with_tolerance <- function(x, y, tolerance = 1e-10) {
  # For each element in x, check if there's at least one element in y that's close enough
  all(sapply(x, function(xi) any(abs(y - xi) < tolerance)))
}

# Test hierarchical beta  
dataTest <- list(rbeta(10, 1, 3), rbeta(10, 1, 3), rbeta(10, 3, 5), rbeta(10, 4, 5), rbeta(10, 6, 3))
dpobjlistTest <- DirichletProcessHierarchicalBeta(dataTest, 1)

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
global_nu <- c(dpobjlistTest$globalParameters[[2]])

cat("\nGlobal parameters:\n")
cat("  mu values:", paste(global_mu, collapse=", "), "\n")
cat("  nu values:", paste(global_nu, collapse=", "), "\n")

# Check each individual DP
for (i in seq_along(dpobjlistTest$indDP)) {
  ind_mu <- c(dpobjlistTest$indDP[[i]]$clusterParameters[[1]])
  ind_nu <- c(dpobjlistTest$indDP[[i]]$clusterParameters[[2]])
  
  cat("\nIndividual DP", i, ":\n")
  cat("  mu values:", paste(ind_mu, collapse=", "), "\n") 
  cat("  nu values:", paste(ind_nu, collapse=", "), "\n")
  
  # Test tolerance function
  mu_result <- all_in_with_tolerance(ind_mu, global_mu, tolerance = 1e-8)
  nu_result <- all_in_with_tolerance(ind_nu, global_nu, tolerance = 1e-8)
  
  cat("  mu tolerance check:", mu_result, "\n")
  cat("  nu tolerance check:", nu_result, "\n")
  
  if (!mu_result) {
    cat("  mu differences:\n")
    for (j in seq_along(ind_mu)) {
      min_diff <- min(abs(global_mu - ind_mu[j]))
      cat("    ind_mu[", j, "] =", ind_mu[j], ", min global diff =", min_diff, "\n")
    }
  }
  
  if (!nu_result) {
    cat("  nu differences:\n")
    for (j in seq_along(ind_nu)) {
      min_diff <- min(abs(global_nu - ind_nu[j]))
      cat("    ind_nu[", j, "] =", ind_nu[j], ", min global diff =", min_diff, "\n")
    }
  }
}