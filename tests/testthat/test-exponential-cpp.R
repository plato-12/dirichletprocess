# tests/testthat/test-exponential-cpp.R
context("Exponential Distribution C++ Implementation")

test_that("Exponential PriorDraw C++ matches R implementation", {
  set.seed(42)
  priorParams <- c(0.01, 0.01)  # alpha0, beta0
  n <- 10

  # Create mixing distribution object for R implementation
  mdObj <- ExponentialMixtureCreate(priorParams)

  # R implementation
  set.seed(42)
  r_result <- PriorDraw(mdObj, n)

  # C++ implementation
  set.seed(42)
  cpp_result <- exponential_prior_draw_cpp(priorParams, n)

  # Check structure
  expect_equal(names(cpp_result), c("lambda"))
  expect_equal(dim(cpp_result$lambda), c(1, 1, n))

  # Check that values are positive (rate parameters must be positive)
  expect_true(all(cpp_result$lambda > 0))

  # Compare statistical properties
  expect_equal(mean(cpp_result$lambda), mean(r_result[[1]]), tolerance = 0.1)
  expect_equal(var(c(cpp_result$lambda)), var(c(r_result[[1]])), tolerance = 0.2)
})

test_that("Exponential PosteriorParameters C++ matches R implementation", {
  priorParams <- c(2, 4)
  x <- matrix(rexp(20, rate = 2), ncol = 1)

  # Create mixing distribution object for R implementation
  mdObj <- ExponentialMixtureCreate(priorParams)

  # R implementation
  r_post_params <- PosteriorParameters(mdObj, x)

  # C++ implementation
  cpp_post_params <- exponential_posterior_parameters_cpp(priorParams, x)

  # Check exact equality (this is deterministic calculation)
  expect_equal(cpp_post_params, r_post_params, tolerance = 1e-10)
})

test_that("Exponential PosteriorDraw C++ matches R implementation", {
  set.seed(123)
  priorParams <- c(2, 4)
  x <- matrix(rexp(25, rate = 3), ncol = 1)
  n <- 5

  # Create mixing distribution object for R implementation
  mdObj <- ExponentialMixtureCreate(priorParams)

  # R implementation
  set.seed(123)
  r_result <- PosteriorDraw(mdObj, x, n)

  # C++ implementation
  set.seed(123)
  cpp_result <- exponential_posterior_draw_cpp(priorParams, x, n)

  # Check structure
  expect_equal(names(cpp_result), c("lambda"))
  expect_equal(dim(cpp_result$lambda), c(1, 1, n))

  # Check that values are positive
  expect_true(all(cpp_result$lambda > 0))

  # Since we use the same seed, results should be very close
  expect_equal(mean(cpp_result$lambda), mean(r_result[[1]]), tolerance = 0.05)
})

test_that("Exponential Likelihood C++ matches R implementation", {
  x <- seq(0.1, 2, by = 0.2)
  lambda <- 2.5

  # Create theta list as expected by R implementation
  theta <- list(array(lambda, dim = c(1, 1, 1)))

  # R implementation
  mdObj <- ExponentialMixtureCreate(c(0.01, 0.01))
  r_lik <- Likelihood(mdObj, x, theta)

  # C++ implementation
  cpp_lik <- exponential_likelihood_cpp(x, lambda)

  # Should be exactly equal
  expect_equal(cpp_lik, r_lik, tolerance = 1e-10)
})

test_that("Exponential Predictive C++ matches R implementation", {
  set.seed(456)
  priorParams <- c(2, 4)
  x <- rexp(5, rate = 2)

  # Create mixing distribution object for R implementation
  mdObj <- ExponentialMixtureCreate(priorParams)

  # R implementation
  r_pred <- Predictive(mdObj, x)

  # C++ implementation
  cpp_pred <- exponential_predictive_cpp(priorParams, x)

  # Should be exactly equal
  expect_equal(cpp_pred, r_pred, tolerance = 1e-10)
})

test_that("ConjugateExponentialClusterComponentUpdate C++ works correctly", {
  set.seed(789)
  # Create a simple DP object with exponential data
  y <- c(rexp(15, rate = 5), rexp(15, rate = 1))
  dp <- DirichletProcessExponential(y)

  # Manually set up clusters
  dp$clusterLabels <- c(rep(1, 15), rep(2, 15))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(15, 15)
  dp$clusterParameters <- list(
    array(c(5, 1), dim = c(1, 1, 2))  # Initial rate parameters
  )

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run C++ update
  set.seed(999)
  cpp_result <- conjugate_exponential_cluster_component_update_cpp(dp_cpp)

  # Convert back to 1-indexed
  cpp_result$clusterLabels <- cpp_result$clusterLabels + 1

  # Check that the structure is preserved
  expect_true(all(cpp_result$clusterLabels %in% 1:30))
  expect_equal(sum(cpp_result$pointsPerCluster), 30)
  expect_true(cpp_result$numberClusters >= 1)
  expect_true(cpp_result$numberClusters <= 30)
})

test_that("ConjugateExponentialClusterParameterUpdate C++ works correctly", {
  set.seed(111)
  # Create data with known structure
  y1 <- rexp(20, rate = 5)
  y2 <- rexp(20, rate = 0.5)
  y <- c(y1, y2)

  dp <- DirichletProcessExponential(y)

  # Set up known clustering
  dp$clusterLabels <- c(rep(1, 20), rep(2, 20))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(20, 20)
  dp$clusterParameters <- list(
    array(c(1, 1), dim = c(1, 1, 2))  # Initial rate parameters
  )

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run parameter update
  set.seed(222)
  cpp_params <- conjugate_exponential_cluster_parameter_update_cpp(dp_cpp)

  # Check that parameters were updated reasonably
  # Cluster 1 should have rate around 5
  expect_true(abs(cpp_params[[1]][1] - 5) < 2)
  # Cluster 2 should have rate around 0.5
  expect_true(abs(cpp_params[[1]][2] - 0.5) < 0.5)

  # All rates should be positive
  expect_true(all(cpp_params[[1]] > 0))
})

test_that("End-to-end Exponential C++ sampler test", {
  set.seed(2025)
  # Generate data with clear clusters
  y <- c(
    rexp(25, rate = 10),   # Fast decay
    rexp(25, rate = 2),    # Medium decay
    rexp(25, rate = 0.5)   # Slow decay
  )

  # Shuffle the data
  y <- sample(y)

  # Initialize DP
  dp <- DirichletProcessExponential(y)

  # Convert to format expected by C++ (0-indexed clusters)
  dp$clusterLabels <- dp$clusterLabels - 1

  # Store initial state
  initial_clusters <- length(unique(dp$clusterLabels))

  # Run 20 iterations of C++ sampler
  cluster_history <- numeric(20)

  for (iter in 1:20) {
    # Update cluster assignments
    update_result <- conjugate_exponential_cluster_component_update_cpp(dp)

    dp$clusterLabels <- update_result$clusterLabels
    dp$pointsPerCluster <- update_result$pointsPerCluster
    dp$numberClusters <- update_result$numberClusters
    dp$clusterParameters <- update_result$clusterParameters

    # Update cluster parameters
    dp$clusterParameters <- conjugate_exponential_cluster_parameter_update_cpp(dp)

    # Record number of clusters
    cluster_history[iter] <- dp$numberClusters

    # Basic sanity checks after each iteration
    expect_equal(length(dp$clusterLabels), length(y))
    expect_equal(sum(dp$pointsPerCluster), length(y))
    expect_true(dp$numberClusters >= 1)
    expect_true(dp$numberClusters <= length(y))
    expect_equal(length(dp$clusterParameters[[1]]), dp$numberClusters)
  }

  # Check that we found reasonable clusters (should be around 3)
  final_clusters <- dp$numberClusters
  expect_true(final_clusters >= 2 && final_clusters <= 5)

  # Convert back to 1-indexed for inspection
  final_labels <- dp$clusterLabels + 1

  # Print summary for manual inspection
  cat("\nEnd-to-end Exponential test summary:\n")
  cat("Initial clusters:", initial_clusters, "\n")
  cat("Final clusters:", final_clusters, "\n")
  cat("Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
  cat("Cluster rates:", round(as.numeric(dp$clusterParameters[[1]]), 2), "\n")

  # Verify that the rates make sense (sorted)
  rates <- sort(as.numeric(dp$clusterParameters[[1]]))
  if (final_clusters == 3) {
    # If we found 3 clusters, they should roughly correspond to our true rates
    expect_true(rates[1] < 1.5)    # Slowest rate
    expect_true(rates[2] > 1 && rates[2] < 5)  # Medium rate
    expect_true(rates[3] > 5)       # Fastest rate
  }
})

test_that("Performance comparison: Exponential R vs C++", {
  skip_if_not(interactive(), "Performance test only run interactively")

  set.seed(3000)
  y <- rexp(100, rate = 2)

  # R implementation
  dp_r <- DirichletProcessExponential(y)
  time_r <- system.time({
    for (i in 1:10) {
      dp_r <- ClusterComponentUpdate(dp_r)
      dp_r <- ClusterParameterUpdate(dp_r)
    }
  })

  # C++ implementation
  dp_cpp <- DirichletProcessExponential(y)
  dp_cpp$clusterLabels <- dp_cpp$clusterLabels - 1

  time_cpp <- system.time({
    for (i in 1:10) {
      update_result <- conjugate_exponential_cluster_component_update_cpp(dp_cpp)
      dp_cpp$clusterLabels <- update_result$clusterLabels
      dp_cpp$pointsPerCluster <- update_result$pointsPerCluster
      dp_cpp$numberClusters <- update_result$numberClusters
      dp_cpp$clusterParameters <- update_result$clusterParameters

      dp_cpp$clusterParameters <- conjugate_exponential_cluster_parameter_update_cpp(dp_cpp)
    }
  })

  cat("\nExponential Performance Comparison (10 iterations):\n")
  cat("R implementation:", round(time_r["elapsed"], 3), "seconds\n")
  cat("C++ implementation:", round(time_cpp["elapsed"], 3), "seconds\n")
  cat("Speedup:", round(time_r["elapsed"] / time_cpp["elapsed"], 1), "x\n")
})
