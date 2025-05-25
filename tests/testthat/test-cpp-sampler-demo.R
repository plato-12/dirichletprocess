context("C++ Sampler End-to-End Demonstration")

test_that("Complete C++ sampler demonstration", {
  # Generate synthetic data with 3 clear clusters
  set.seed(2025)
  true_means <- c(-5, 0, 5)
  true_sds <- c(0.5, 0.5, 0.5)
  n_per_cluster <- 30

  y <- c(
    rnorm(n_per_cluster, true_means[1], true_sds[1]),
    rnorm(n_per_cluster, true_means[2], true_sds[2]),
    rnorm(n_per_cluster, true_means[3], true_sds[3])
  )

  # Shuffle the data
  y <- sample(y)

  cat("\n=== C++ Sampler Demonstration ===\n")
  cat("Data: 3 clusters with means at", true_means, "\n")
  cat("Total observations:", length(y), "\n\n")

  # Initialize Dirichlet Process
  dp <- DirichletProcessGaussian(y)
  cat("Initial state: 1 cluster\n")

  # Run the C++ sampler for 20 iterations
  n_iterations <- 20
  cluster_history <- numeric(n_iterations)

  # Convert to 0-indexed for C++
  dp$clusterLabels <- dp$clusterLabels - 1

  cat("\nRunning C++ sampler...\n")
  for (iter in 1:n_iterations) {
    # Step 1: Update cluster assignments (Chinese Restaurant Process)
    update_result <- conjugate_cluster_component_update_cpp(dp)

    dp$clusterLabels <- update_result$clusterLabels
    dp$pointsPerCluster <- update_result$pointsPerCluster
    dp$numberClusters <- update_result$numberClusters
    dp$clusterParameters <- update_result$clusterParameters

    # Step 2: Update cluster parameters
    dp$clusterParameters <- conjugate_cluster_parameter_update_cpp(dp)

    # Record number of clusters
    cluster_history[iter] <- dp$numberClusters

    if (iter %% 5 == 0) {
      cat("  Iteration", iter, ": ", dp$numberClusters, "clusters\n")
    }
  }

  # Final results
  cat("\nFinal results:\n")
  cat("  Number of clusters:", dp$numberClusters, "\n")
  cat("  Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
  cat("  Cluster means:", round(as.numeric(dp$clusterParameters[[1]]), 2), "\n")
  cat("  Cluster sds:", round(as.numeric(dp$clusterParameters[[2]]), 2), "\n")

  # Plot cluster history
  if (requireNamespace("graphics", quietly = TRUE)) {
    plot(1:n_iterations, cluster_history, type = "l",
         xlab = "Iteration", ylab = "Number of Clusters",
         main = "Cluster Evolution (C++ Sampler)",
         ylim = c(0, max(cluster_history) + 1))
    abline(h = 3, col = "red", lty = 2)
    legend("topright", "True clusters = 3", col = "red", lty = 2)
  }

  # Validate results
  expect_true(dp$numberClusters >= 2 && dp$numberClusters <= 5,
              "Should find approximately 3 clusters")
  expect_equal(sum(dp$pointsPerCluster), length(y), info = "All points should be assigned")
  expect_true(all(dp$clusterParameters[[2]] > 0),
              "All standard deviations should be positive")

  # Check that cluster means are reasonable
  cluster_means <- sort(as.numeric(dp$clusterParameters[[1]]))
  if (dp$numberClusters == 3) {
    # If we found 3 clusters, they should be near the true means
    expect_true(abs(cluster_means[1] - true_means[1]) < 1)
    expect_true(abs(cluster_means[2] - true_means[2]) < 1)
    expect_true(abs(cluster_means[3] - true_means[3]) < 1)
  }

  cat("\n=== Demonstration Complete ===\n")
})

test_that("Performance comparison: C++ vs R", {
  skip_if_not(interactive(), "Performance test only run interactively")

  set.seed(3000)
  y <- rnorm(100)

  # R implementation
  dp_r <- DirichletProcessGaussian(y)
  time_r <- system.time({
    for (i in 1:10) {
      dp_r <- ClusterComponentUpdate(dp_r)
      dp_r <- ClusterParameterUpdate(dp_r)
    }
  })

  # C++ implementation
  dp_cpp <- DirichletProcessGaussian(y)
  dp_cpp$clusterLabels <- dp_cpp$clusterLabels - 1

  time_cpp <- system.time({
    for (i in 1:10) {
      update_result <- conjugate_cluster_component_update_cpp(dp_cpp)
      dp_cpp$clusterLabels <- update_result$clusterLabels
      dp_cpp$pointsPerCluster <- update_result$pointsPerCluster
      dp_cpp$numberClusters <- update_result$numberClusters
      dp_cpp$clusterParameters <- update_result$clusterParameters

      dp_cpp$clusterParameters <- conjugate_cluster_parameter_update_cpp(dp_cpp)
    }
  })

  cat("\nPerformance Comparison (10 iterations):\n")
  cat("R implementation:", round(time_r["elapsed"], 3), "seconds\n")
  cat("C++ implementation:", round(time_cpp["elapsed"], 3), "seconds\n")
  cat("Speedup:", round(time_r["elapsed"] / time_cpp["elapsed"], 1), "x\n")
})
