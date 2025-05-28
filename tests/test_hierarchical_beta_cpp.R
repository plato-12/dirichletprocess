# Test script for Hierarchical Beta DP C++ implementation

library(dirichletprocess)

# Set random seed for reproducibility
set.seed(123)

# Generate test data
n_datasets <- 3
n_points_per_dataset <- 20

# Generate data from different Beta distributions
dataList <- list(
  rbeta(n_points_per_dataset, 2, 5),
  rbeta(n_points_per_dataset, 5, 2),
  rbeta(n_points_per_dataset, 3, 3)
)

# Scale data to [0, 1] range
maxY <- 1

# Create hierarchical Beta DP with R implementation
message("Creating hierarchical Beta DP with R implementation...")
dp_r <- DirichletProcessHierarchicalBeta(
  dataList = dataList,
  maxY = maxY,
  priorParameters = c(2, 8),
  hyperPriorParameters = c(1, 0.125),
  gammaPriors = c(2, 4),
  alphaPriors = c(2, 4),
  mhStepSize = c(0.1, 0.1),
  numSticks = 20,
  mhDraws = 10
)

# Fit with R implementation
message("Fitting with R implementation...")
system.time(dp_r_fit <- Fit(dp_r, its = 10, progressBar = FALSE))

# Enable C++ implementations
enable_cpp_hierarchical_samplers(TRUE)

# Create hierarchical Beta DP with C++ implementation
message("\nCreating hierarchical Beta DP with C++ implementation...")
dp_cpp <- DirichletProcessHierarchicalBeta(
  dataList = dataList,
  maxY = maxY,
  priorParameters = c(2, 8),
  hyperPriorParameters = c(1, 0.125),
  gammaPriors = c(2, 4),
  alphaPriors = c(2, 4),
  mhStepSize = c(0.1, 0.1),
  numSticks = 20,
  mhDraws = 10
)

# Fit with C++ implementation
message("Fitting with C++ implementation...")
system.time(dp_cpp_fit <- Fit(dp_cpp, its = 10, progressBar = FALSE))

# Compare results
message("\nComparing results...")
message(sprintf("R implementation - Final gamma: %.4f", dp_r_fit$gamma))
message(sprintf("C++ implementation - Final gamma: %.4f", dp_cpp_fit$gamma))

message(sprintf("\nR implementation - Number of global parameters: %d",
                length(dp_r_fit$globalParameters[[1]])))
message(sprintf("C++ implementation - Number of global parameters: %d",
                length(dp_cpp_fit$globalParameters[[1]])))

# Compare individual DP properties
for (i in seq_along(dataList)) {
  message(sprintf("\nDataset %d:", i))
  message(sprintf("  R - Number of clusters: %d, alpha: %.4f",
                  dp_r_fit$indDP[[i]]$numberClusters,
                  dp_r_fit$indDP[[i]]$alpha))
  message(sprintf("  C++ - Number of clusters: %d, alpha: %.4f",
                  dp_cpp_fit$indDP[[i]]$numberClusters,
                  dp_cpp_fit$indDP[[i]]$alpha))
}

# Test individual update functions
message("\nTesting individual update functions...")

# Test ClusterComponentUpdate
dp_test <- dp_cpp
dp_test_updated <- ClusterComponentUpdate(dp_test)
message("ClusterComponentUpdate: OK")

# Test GlobalParameterUpdate
dp_test_updated <- GlobalParameterUpdate(dp_test)
message("GlobalParameterUpdate: OK")

# Test UpdateG0
dp_test_updated <- UpdateG0(dp_test)
message("UpdateG0: OK")

# Test UpdateGamma
dp_test_updated <- UpdateGamma(dp_test)
message("UpdateGamma: OK")

message("\nAll tests completed successfully!")

# Disable C++ implementations
enable_cpp_hierarchical_samplers(FALSE)
