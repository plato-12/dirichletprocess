# Debug hierarchical beta aux parameter structure
library(dirichletprocess)

# Create simple test data
group1_data <- rbeta(10, 2, 5)
group2_data <- rbeta(10, 5, 2)
hierarchical_data <- list(group1_data, group2_data)

# Create hierarchical beta DP
dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)

# Check structure of individual DPs
cat("Individual DP 1 class:", class(dp$indDP[[1]]), "\n")
cat("Individual DP 1 mixing distribution class:", class(dp$indDP[[1]]$mixingDistribution), "\n")

# Check if aux exists and its structure
if ("aux" %in% names(dp$indDP[[1]])) {
  cat("DP 1 has aux parameters\n")
  cat("aux length:", length(dp$indDP[[1]]$aux), "\n")
  cat("aux[[1]] names:", names(dp$indDP[[1]]$aux[[1]]), "\n")
  cat("aux[[1]]$mu structure:", str(dp$indDP[[1]]$aux[[1]]$mu), "\n")
  cat("aux[[1]]$nu structure:", str(dp$indDP[[1]]$aux[[1]]$nu), "\n")
} else {
  cat("DP 1 does NOT have aux parameters\n")
}

# Check cluster parameters structure
cat("Cluster parameters names:", names(dp$indDP[[1]]$clusterParameters), "\n")
cat("Cluster parameters mu structure:", str(dp$indDP[[1]]$clusterParameters$mu), "\n")
cat("Cluster parameters nu structure:", str(dp$indDP[[1]]$clusterParameters$nu), "\n")