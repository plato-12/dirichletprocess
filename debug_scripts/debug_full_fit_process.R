# Debug full hierarchical beta fit process
library(dirichletprocess)

# Disable C++ to isolate R implementation issue
set_use_cpp(FALSE)

# Create test data similar to the failing test
group1_data <- rbeta(25, 2, 5)
group2_data <- rbeta(25, 5, 2) 
group3_data <- rbeta(25, 1, 1)
hierarchical_data <- list(group1_data, group2_data, group3_data)

# Create hierarchical beta DP
dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)

cat("Starting manual MCMC steps...\n")

for (i in 1:3) {
  cat("Iteration", i, "\n")
  
  # Step 1: ClusterComponentUpdate
  tryCatch({
    dp <- ClusterComponentUpdate(dp)
    cat("  ClusterComponentUpdate succeeded\n")
  }, error = function(e) {
    cat("  ClusterComponentUpdate failed:", e$message, "\n")
    return()
  })
  
  # Step 2: ClusterParameterUpdate  
  tryCatch({
    dp <- ClusterParameterUpdate(dp)
    cat("  ClusterParameterUpdate succeeded\n")
  }, error = function(e) {
    cat("  ClusterParameterUpdate failed:", e$message, "\n")
    return()
  })
  
  # Step 3: UpdateAlpha
  tryCatch({
    dp <- UpdateAlpha(dp)
    cat("  UpdateAlpha succeeded\n")
  }, error = function(e) {
    cat("  UpdateAlpha failed:", e$message, "\n")
    return()
  })
  
  # Step 4: GlobalParameterUpdate
  tryCatch({
    dp <- GlobalParameterUpdate(dp)
    cat("  GlobalParameterUpdate succeeded\n")
  }, error = function(e) {
    cat("  GlobalParameterUpdate failed:", e$message, "\n")
    
    # Debug when it fails
    cat("  Debug info for iteration", i, ":\n")
    cat("  Number of clusters in DP 1:", dp$indDP[[1]]$numberClusters, "\n")
    cat("  Number of clusters in DP 2:", dp$indDP[[2]]$numberClusters, "\n") 
    cat("  Number of clusters in DP 3:", dp$indDP[[3]]$numberClusters, "\n")
    
    # Check parameter dimensions
    for (j in 1:3) {
      cat("  DP", j, "clusterParameters[[1]] dims:", dim(dp$indDP[[j]]$clusterParameters[[1]]), "\n")
      cat("  DP", j, "clusterParameters[[2]] dims:", dim(dp$indDP[[j]]$clusterParameters[[2]]), "\n")
    }
    
    return()
  })
}

cat("Manual MCMC complete\n")