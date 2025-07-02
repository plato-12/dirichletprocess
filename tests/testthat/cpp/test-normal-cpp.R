context("Normal Distribution C++ Implementation")

test_that("Normal PriorDraw C++ matches R implementation", {
  set.seed(42)
  priorParams <- c(0, 1, 1, 1)  # mu0, kappa0, alpha0, beta0
  n <- 10

  # Create mixing distribution object for R implementation
  mdObj <- MixingDistribution("normal", priorParams, "conjugate")

  # R implementation
  set.seed(42)
  r_result <- PriorDraw(mdObj, n)

  # C++ implementation
  set.seed(42)
  cpp_result <- normal_prior_draw_cpp(priorParams, n)

  # Check structure
  expect_equal(names(cpp_result), names(r_result))
  expect_equal(dim(cpp_result[[1]]), dim(r_result[[1]]))
  expect_equal(dim(cpp_result[[2]]), dim(r_result[[2]]))

  # Check values are statistically similar (not exact due to RNG differences)
  # Compare means
  expect_equal(mean(cpp_result[[1]]), mean(r_result[[1]]), tolerance = 0.1)
  expect_equal(mean(cpp_result[[2]]), mean(r_result[[2]]), tolerance = 0.1)

  # Compare variances
  expect_equal(var(c(cpp_result[[1]])), var(c(r_result[[1]])), tolerance = 0.2)
  expect_equal(var(c(cpp_result[[2]])), var(c(r_result[[2]])), tolerance = 0.2)
})

test_that("Normal PosteriorDraw C++ matches R implementation", {
  set.seed(123)
  priorParams <- c(0, 1, 1, 1)
  x <- matrix(rnorm(20), ncol = 1)
  n <- 5

  # Create mixing distribution object for R implementation
  mdObj <- MixingDistribution("normal", priorParams, "conjugate")

  # R implementation
  set.seed(123)
  r_result <- PosteriorDraw(mdObj, x, n)

  # C++ implementation
  set.seed(123)
  cpp_result <- normal_posterior_draw_cpp(priorParams, x, n)

  # Check structure
  expect_equal(names(cpp_result), names(r_result))
  expect_equal(dim(cpp_result[[1]]), dim(r_result[[1]]))
  expect_equal(dim(cpp_result[[2]]), dim(r_result[[2]]))

  # Check posterior parameters are calculated correctly
  r_post_params <- PosteriorParameters(mdObj, x)
  cpp_post_params <- normal_posterior_parameters_cpp(priorParams, x)

  expect_equal(cpp_post_params, r_post_params, tolerance = 1e-10)
})

test_that("ConjugateClusterComponentUpdate C++ matches behavior", {
  set.seed(456)
  # Create a simple DP object
  y <- c(rnorm(10, -2, 0.5), rnorm(10, 2, 0.5))
  dp <- DirichletProcessGaussian(y)

  # Manually set up clusters
  dp$clusterLabels <- c(rep(1, 10), rep(2, 10))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(10, 10)
  dp$clusterParameters <- list(
    array(c(-2, 2), dim = c(1, 1, 2)),
    array(c(0.5, 0.5), dim = c(1, 1, 2))
  )

  # Convert cluster labels to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run C++ update
  set.seed(789)
  cpp_result <- conjugate_cluster_component_update_cpp(dp_cpp)

  # Convert back to 1-indexed
  cpp_result$clusterLabels <- cpp_result$clusterLabels + 1

  # Check that the structure is preserved
  expect_true(all(cpp_result$clusterLabels %in% 1:20))
  expect_equal(sum(cpp_result$pointsPerCluster), 20)
  expect_true(cpp_result$numberClusters >= 1)
  expect_true(cpp_result$numberClusters <= 20)
})

test_that("ConjugateClusterParameterUpdate C++ works correctly", {
  set.seed(789)
  # Create a simple DP object with known clusters
  y1 <- rnorm(15, mean = -3, sd = 0.5)
  y2 <- rnorm(15, mean = 3, sd = 0.5)
  y <- c(y1, y2)

  dp <- DirichletProcessGaussian(y)

  # Set up known clustering
  dp$clusterLabels <- c(rep(1, 15), rep(2, 15))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(15, 15)
  dp$clusterParameters <- list(
    array(c(0, 0), dim = c(1, 1, 2)),  # Initial means
    array(c(1, 1), dim = c(1, 1, 2))   # Initial sds
  )

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run parameter update
  set.seed(999)
  cpp_params <- conjugate_cluster_parameter_update_cpp(dp_cpp)

  # Check that parameters were updated reasonably
  # Cluster 1 should have mean around -3
  expect_true(abs(cpp_params[[1]][1] - (-3)) < 1)
  # Cluster 2 should have mean around 3
  expect_true(abs(cpp_params[[1]][2] - 3) < 1)

  # Standard deviations should be reasonable
  expect_true(all(cpp_params[[2]] > 0))
  expect_true(all(cpp_params[[2]] < 2))
})

test_that("End-to-end C++ sampler test", {
  set.seed(1234)
  # Generate data with clear clusters
  y <- c(rnorm(20, -3, 0.5), rnorm(20, 0, 0.5), rnorm(20, 3, 0.5))

  # Initialize DP
  dp <- DirichletProcessGaussian(y)

  # Convert to format expected by C++ (0-indexed clusters)
  dp$clusterLabels <- dp$clusterLabels - 1

  # Store initial state
  initial_clusters <- length(unique(dp$clusterLabels))

  # Run 10 iterations of C++ sampler
  for (iter in 1:10) {
    # Update cluster assignments
    update_result <- conjugate_cluster_component_update_cpp(dp)

    # Update dp object with results
    dp$clusterLabels <- update_result$clusterLabels
    dp$pointsPerCluster <- update_result$pointsPerCluster
    dp$numberClusters <- update_result$numberClusters
    dp$clusterParameters <- update_result$clusterParameters

    # Update cluster parameters
    dp$clusterParameters <- conjugate_cluster_parameter_update_cpp(dp)

    # Basic sanity checks after each iteration
    expect_equal(length(dp$clusterLabels), length(y))
    expect_equal(sum(dp$pointsPerCluster), length(y))
    expect_true(dp$numberClusters >= 1)
    expect_true(dp$numberClusters <= length(y))
    expect_equal(length(dp$clusterParameters[[1]]), dp$numberClusters)
    expect_equal(length(dp$clusterParameters[[2]]), dp$numberClusters)
  }

  # Check that we found reasonable clusters (should be around 3)
  final_clusters <- dp$numberClusters
  expect_true(final_clusters >= 2 && final_clusters <= 8)

  # Convert back to 1-indexed for inspection
  final_labels <- dp$clusterLabels + 1

  # Print summary for manual inspection
  cat("\nEnd-to-end test summary:\n")
  cat("Initial clusters:", initial_clusters, "\n")
  cat("Final clusters:", final_clusters, "\n")
  cat("Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
  cat("Cluster means:", as.numeric(dp$clusterParameters[[1]]), "\n")
})
