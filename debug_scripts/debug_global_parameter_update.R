# Debug GlobalParameterUpdate for hierarchical beta
library(dirichletprocess)

# Disable C++ to isolate R implementation issue
set_use_cpp(FALSE)

# Create simple test data
group1_data <- rbeta(5, 2, 5)
group2_data <- rbeta(5, 5, 2)
hierarchical_data <- list(group1_data, group2_data)

# Create hierarchical beta DP
dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)

cat("Testing GlobalParameterUpdate...\n")
cat("Individual DP 1 clusterParameters structure:\n")
str(dp$indDP[[1]]$clusterParameters)

tryCatch({
  dp_updated <- GlobalParameterUpdate(dp)
  cat("GlobalParameterUpdate succeeded\n")
}, error = function(e) {
  cat("GlobalParameterUpdate failed with error:", e$message, "\n")
  
  # Debug the parameter structures
  cat("\nDiagnostics:\n")
  cat("Global parameters structure:\n")
  str(dp$globalParameters)
  
  cat("\nIndividual DP 1 cluster parameters structure:\n")
  str(dp$indDP[[1]]$clusterParameters)
  
  # Check dimensions
  cat("\nDimensions:\n")
  cat("globalParameters[[1]] dims:", dim(dp$globalParameters[[1]]), "\n")
  cat("globalParameters[[2]] dims:", dim(dp$globalParameters[[2]]), "\n")
  cat("indDP[[1]]$clusterParameters[[1]] dims:", dim(dp$indDP[[1]]$clusterParameters[[1]]), "\n")
  cat("indDP[[1]]$clusterParameters[[2]] dims:", dim(dp$indDP[[1]]$clusterParameters[[2]]), "\n")
})

cat("Debugging complete\n")