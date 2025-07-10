# test_beta_fix.R
library(dirichletprocess)

# Test function with detailed output
test_beta_dp <- function(n = 20, iterations = 10, seed = 123) {
  set.seed(seed)

  # Generate test data
  data <- rbeta(n, 2, 5)

  # Create and initialize DP
  dp <- DirichletProcessBeta(data, verbose = FALSE)

  cat("Initial state:\n")
  cat("  Number of clusters:", dp$numberClusters, "\n")
  cat("  Points per cluster:", dp$pointsPerCluster, "\n")
  cat("  Sum points:", sum(dp$pointsPerCluster), "\n")
  cat("  Expected:", n, "\n\n")

  # Run iterations
  for (iter in 1:iterations) {
    dp <- ClusterComponentUpdate(dp)

    # Validate state
    sum_points <- sum(dp$pointsPerCluster)
    if (sum_points != n) {
      cat("ERROR at iteration", iter, ":\n")
      cat("  Sum points:", sum_points, "Expected:", n, "\n")
      cat("  Clusters:", dp$numberClusters, "\n")
      cat("  Points per cluster:", dp$pointsPerCluster, "\n")
      cat("  Label counts:", table(dp$clusterLabels), "\n")
      stop("Point count mismatch!")
    }

    if (any(dp$clusterLabels <= 0)) {
      stop("Non-positive cluster labels found!")
    }

    if (iter %% 5 == 0) {
      cat("Iteration", iter, "- OK (clusters:", dp$numberClusters, ")\n")
    }
  }

  cat("\nFinal state:\n")
  cat("  Number of clusters:", dp$numberClusters, "\n")
  cat("  Points per cluster:", dp$pointsPerCluster, "\n")
  cat("  All cluster labels > 0:", all(dp$clusterLabels > 0), "\n")
  cat("  Sum points correct:", sum(dp$pointsPerCluster) == n, "\n")

  return(dp)
}

# Run the test
result <- test_beta_dp(n = 50, iterations = 20)
