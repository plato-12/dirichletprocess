# tests/testthat/test-exponential-cpp.R
context("Exponential Distribution C++ Implementation")

test_that("Exponential PriorDraw C++ matches R implementation", {
  set.seed(42)
  priorParams <- c(0.01, 0.01)  # alpha0, beta0
  n <- 10

  # R implementation - implement directly to avoid S3 dispatch issues
  set.seed(42)
  r_draws <- rgamma(n, priorParams[1], priorParams[2])
  r_result <- list(array(r_draws, dim=c(1,1,n)))

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

  # R implementation - implement directly
  alpha_n <- priorParams[1] + length(x)
  beta_n <- priorParams[2] + sum(x)
  r_post_params <- matrix(c(alpha_n, beta_n), nrow = 1)

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

  # R implementation - implement directly
  set.seed(123)
  r_theta <- rgamma(n, priorParams[1] + length(x), priorParams[2] + sum(x))
  r_result <- list(array(r_theta, dim=c(1,1,n)))

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

  # R implementation - implement directly
  r_lik <- dexp(x, lambda)

  # C++ implementation
  cpp_lik <- exponential_likelihood_cpp(x, lambda)

  # Should be exactly equal
  expect_equal(cpp_lik, r_lik, tolerance = 1e-10)
})

test_that("Exponential Predictive C++ matches R implementation", {
  set.seed(456)
  priorParams <- c(2, 4)
  x <- rexp(5, rate = 2)

  # R implementation - implement directly based on exponential_gamma.R
  r_pred <- numeric(length(x))
  for(i in seq_along(x)){
    alphaPost <- priorParams[1] + 1  # length(x[i]) = 1
    betaPost <- priorParams[2] + x[i]  # sum(x[i]) = x[i]
    r_pred[i] <- (gamma(alphaPost)/gamma(priorParams[1])) *
      ((priorParams[2]^priorParams[1])/((betaPost)^alphaPost))
  }

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

  # Initialize DP properly
  dp <- DirichletProcessExponential(y, alphaPriors = c(2, 2))

  # Set higher initial alpha to encourage more clusters
  dp$alpha <- 2.0

  # Ensure predictiveArray is initialized if not already
  if (is.null(dp$predictiveArray) || length(dp$predictiveArray) == 0) {
    dp$predictiveArray <- Predictive(dp$mixingDistribution, dp$data)
  }

  # Convert to format expected by C++ (0-indexed clusters)
  dp$clusterLabels <- dp$clusterLabels - 1

  # Also ensure alphaPriorParameters is set for C++ updateAlpha
  dp$alphaPriorParameters <- c(2, 2)

  # Store initial state
  initial_clusters <- length(unique(dp$clusterLabels))

  # Run 50 iterations of C++ sampler (more iterations for better mixing)
  cluster_history <- numeric(50)

  for (iter in 1:50) {
    # Update cluster assignments
    update_result <- conjugate_exponential_cluster_component_update_cpp(dp)

    dp$clusterLabels <- update_result$clusterLabels
    dp$pointsPerCluster <- update_result$pointsPerCluster
    dp$numberClusters <- update_result$numberClusters
    dp$clusterParameters <- update_result$clusterParameters

    # Update cluster parameters
    dp$clusterParameters <- conjugate_exponential_cluster_parameter_update_cpp(dp)

    # CRITICAL: Update alpha (concentration parameter)
    # The C++ object needs to be reconstructed with updated state
    dp_for_alpha <- dp
    dp_for_alpha$n <- length(y)

    # Create a new C++ object and update alpha
    # Note: This is a workaround since we can't directly call updateAlpha on the C++ object
    # In a full implementation, you'd have a dedicated C++ function for this
    old_alpha <- dp$alpha
    x <- rbeta(1, dp$alpha + 1, dp$n)
    pi1 <- dp$alphaPriorParameters[1] + dp$numberClusters - 1
    pi2 <- dp$n * (dp$alphaPriorParameters[2] - log(x))
    pi_ratio <- pi1 / (pi1 + pi2)

    if (runif(1) < pi_ratio) {
      postShape <- dp$alphaPriorParameters[1] + dp$numberClusters
    } else {
      postShape <- dp$alphaPriorParameters[1] + dp$numberClusters - 1
    }
    postRate <- dp$alphaPriorParameters[2] - log(x)

    dp$alpha <- rgamma(1, postShape, 1/postRate)

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

  # With proper alpha updates, we should find at least 2 clusters
  expect_true(final_clusters >= 2)
  # But not too many
  expect_true(final_clusters <= 10)

  # Convert back to 1-indexed for inspection
  final_labels <- dp$clusterLabels + 1

  # Print summary for manual inspection
  cat("\nEnd-to-end Exponential test summary:\n")
  cat("Initial clusters:", initial_clusters, "\n")
  cat("Final clusters:", final_clusters, "\n")
  cat("Final alpha:", round(dp$alpha, 2), "\n")
  cat("Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
  cat("Cluster rates:", round(as.numeric(dp$clusterParameters[[1]]), 2), "\n")
  cat("Mean clusters over iterations:", round(mean(cluster_history), 2), "\n")

  # Verify that the rates are reasonable
  rates <- sort(as.numeric(dp$clusterParameters[[1]]))

  # At minimum, we should have separated the slow decay (rate ~0.5) from others
  if (final_clusters >= 2) {
    expect_true(min(rates) < 1.5)  # Should have found the slow decay cluster
    expect_true(max(rates) > 1.5)  # Should have found faster decay cluster(s)
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
      dp_r <- UpdateAlpha(dp_r)  # Include alpha update
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

      # Update alpha
      x <- rbeta(1, dp_cpp$alpha + 1, dp_cpp$n)
      pi1 <- dp_cpp$alphaPriorParameters[1] + dp_cpp$numberClusters - 1
      pi2 <- dp_cpp$n * (dp_cpp$alphaPriorParameters[2] - log(x))
      pi_ratio <- pi1 / (pi1 + pi2)

      if (runif(1) < pi_ratio) {
        postShape <- dp_cpp$alphaPriorParameters[1] + dp_cpp$numberClusters
      } else {
        postShape <- dp_cpp$alphaPriorParameters[1] + dp_cpp$numberClusters - 1
      }
      postRate <- dp_cpp$alphaPriorParameters[2] - log(x)

      dp_cpp$alpha <- rgamma(1, postShape, 1/postRate)
    }
  })

  cat("\nExponential Performance Comparison (10 iterations):\n")
  cat("R implementation:", round(time_r["elapsed"], 3), "seconds\n")
  cat("C++ implementation:", round(time_cpp["elapsed"], 3), "seconds\n")
  cat("Speedup:", round(time_r["elapsed"] / time_cpp["elapsed"], 1), "x\n")
})
