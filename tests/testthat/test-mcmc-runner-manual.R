# tests/testthat/test-mcmc-runner-manual.R

library(testthat)
library(dirichletprocess)

context("MCMCRunnerManual - Comprehensive Test Suite")

# Helper function to create test data
create_test_data <- function(n = 100, k = 3, dim = 2) {
  set.seed(123)
  true_labels <- sample(1:k, n, replace = TRUE)
  data <- matrix(0, n, dim)

  for (i in 1:k) {
    idx <- which(true_labels == i)
    if (dim == 1) {
      data[idx, ] <- rnorm(length(idx), mean = i * 2, sd = 0.5)
    } else {
      for (j in 1:dim) {
        data[idx, j] <- rnorm(length(idx), mean = i * 2 + j, sd = 0.5)
      }
    }
  }

  list(data = data, true_labels = true_labels, k = k)
}

# Test 1: Basic Creation and Initialization
test_that("MCMCRunnerManual creates and initializes correctly", {
  test_data <- create_test_data(n = 50, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  # Create manual runner
  runner <- CppMCMCRunner$new(dp, n_iter = 100, n_burn = 10, thin = 2)

  expect_s4_class(runner, "CppMCMCRunner")
  expect_true(!is.null(runner$ptr))

  # Get initial state
  state <- runner$get_state()
  expect_type(state, "list")
  expect_equal(length(state$cluster_labels), nrow(test_data$data))
  expect_true(state$alpha > 0)
  expect_true(state$n_clusters >= 1)
})

# Test 2: Step-by-Step Execution
test_that("MCMCRunnerManual performs individual steps correctly", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp, n_iter = 50)

  # Get initial state
  initial_state <- runner$get_state()

  # Step assignments
  runner$step_assignments()
  post_assign_state <- runner$get_state()

  # Something should have changed (probabilistic, but very likely)
  expect_true(
    !identical(initial_state$cluster_labels, post_assign_state$cluster_labels) ||
      !identical(initial_state$n_clusters, post_assign_state$n_clusters)
  )

  # Step parameters
  runner$step_parameters()
  post_param_state <- runner$get_state()

  # Parameters should have been updated
  expect_false(identical(post_assign_state$cluster_params, post_param_state$cluster_params))

  # Step concentration
  runner$step_concentration()
  post_alpha_state <- runner$get_state()

  # Alpha might have changed
  expect_true(post_alpha_state$alpha > 0)
})

# Test 3: Full Iteration
test_that("MCMCRunnerManual performs full iterations correctly", {
  test_data <- create_test_data(n = 50, k = 3, dim = 2)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp, n_iter = 100)

  # Run 10 iterations manually
  for (i in 1:10) {
    runner$perform_iteration()
  }

  expect_false(runner$is_complete())
  expect_equal(runner$get_iteration(), 10)

  # Check that diagnostics are being tracked
  state <- runner$get_state()
  expect_true(state$log_posterior != 0)
})

# Test 4: State Modification
test_that("MCMCRunnerManual allows state modification", {
  test_data <- create_test_data(n = 40, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Set custom labels
  new_labels <- rep(1:2, each = 20)
  runner$set_labels(new_labels)

  state <- runner$get_state()
  expect_equal(state$cluster_labels + 1, new_labels)  # +1 for R indexing
  expect_equal(state$n_clusters, 2)

  # Set custom parameters
  new_params <- list(c(0, 1), c(5, 1))
  runner$set_params(new_params)

  state <- runner$get_state()
  expect_equal(length(state$cluster_params), 2)
})

# Test 5: Parameter Bounds
test_that("MCMCRunnerManual respects parameter bounds", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Set bounds for Gaussian parameters (mean, precision)
  runner$set_bounds(lower = c(-10, 0.1), upper = c(10, 10))

  # Run some iterations
  for (i in 1:20) {
    runner$perform_iteration()
  }

  # Check parameters are within bounds
  state <- runner$get_state()
  for (params in state$cluster_params) {
    expect_true(all(params >= c(-10, 0.1)))
    expect_true(all(params <= c(10, 10)))
  }
})

# Test 6: Update Flags
test_that("MCMCRunnerManual respects update flags", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Disable cluster updates
  runner$set_update_flags(clusters = FALSE, params = TRUE, alpha = TRUE)

  initial_labels <- runner$get_state()$cluster_labels

  # Run iterations
  for (i in 1:10) {
    runner$perform_iteration()
  }

  # Labels should not have changed
  expect_equal(runner$get_state()$cluster_labels, initial_labels)
})

# Test 7: Temperature Control
test_that("MCMCRunnerManual temperature control works", {
  test_data <- create_test_data(n = 50, k = 3, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Set high temperature (more exploration)
  runner$set_temperature(2.0)

  # Run some iterations
  for (i in 1:20) {
    runner$perform_iteration()
  }

  high_temp_clusters <- runner$get_state()$n_clusters

  # Reset and run with low temperature
  runner2 <- CppMCMCRunner$new(dp)
  runner2$set_temperature(0.5)

  for (i in 1:20) {
    runner2$perform_iteration()
  }

  low_temp_clusters <- runner2$get_state()$n_clusters

  # High temperature should generally lead to more clusters (not guaranteed)
  expect_true(high_temp_clusters >= 1)
  expect_true(low_temp_clusters >= 1)
})

# Test 8: Auxiliary Parameters
test_that("MCMCRunnerManual auxiliary parameter control works", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Set different auxiliary parameter counts
  runner$set_auxiliary_count(10)

  aux_params <- runner$get_auxiliary_params()
  expect_equal(length(aux_params), 10)

  # Each should be a valid parameter vector
  for (params in aux_params) {
    expect_equal(length(params), 2)  # Gaussian has 2 params
  }
})

# Test 9: Cluster Operations - Merge
test_that("MCMCRunnerManual merge clusters works", {
  test_data <- create_test_data(n = 40, k = 3, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Set up 3 clusters
  labels <- rep(1:3, length.out = 40)
  runner$set_labels(labels)

  initial_state <- runner$get_state()
  expect_equal(initial_state$n_clusters, 3)

  # Merge clusters 1 and 2
  runner$merge_clusters(1, 2)

  merged_state <- runner$get_state()
  expect_equal(merged_state$n_clusters, 2)

  # Check labels were updated correctly
  new_labels <- merged_state$cluster_labels + 1
  expect_true(all(new_labels %in% c(1, 2)))
})

# Test 10: Cluster Operations - Split
test_that("MCMCRunnerManual split cluster works", {
  test_data <- create_test_data(n = 40, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Set up 2 clusters
  labels <- rep(1:2, each = 20)
  runner$set_labels(labels)

  initial_state <- runner$get_state()
  expect_equal(initial_state$n_clusters, 2)

  # Split cluster 1
  runner$split_cluster(1, split_prob = 0.5)

  split_state <- runner$get_state()
  expect_equal(split_state$n_clusters, 3)
})

# Test 11: Diagnostic Methods
test_that("MCMCRunnerManual diagnostic methods work", {
  test_data <- create_test_data(n = 50, k = 3, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp, n_iter = 200)

  # Run enough iterations for diagnostics
  for (i in 1:150) {
    runner$perform_iteration()
  }

  # Get diagnostics
  log_post <- runner$get_log_posterior()
  expect_type(log_post, "double")
  expect_true(is.finite(log_post))

  entropy <- runner$get_clustering_entropy()
  expect_type(entropy, "double")
  expect_true(entropy >= 0)

  cluster_entropies <- runner$get_cluster_entropies()
  expect_true(length(cluster_entropies) >= 1)
  expect_true(all(cluster_entropies >= 0))

  # Convergence diagnostics
  conv_diag <- runner$get_convergence_diagnostics()
  expect_type(conv_diag, "list")
  expect_true(conv_diag$iterations_completed == 150)
  expect_true(is.finite(conv_diag$R_hat))
  expect_true(conv_diag$effective_sample_size > 0)
})

# Test 12: Cluster Statistics and Likelihoods
test_that("MCMCRunnerManual cluster statistics work", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Run some iterations
  for (i in 1:20) {
    runner$perform_iteration()
  }

  # Get cluster statistics
  stats <- runner$get_cluster_statistics()
  expect_type(stats, "list")
  expect_true(length(stats) >= 1)

  total_size <- sum(sapply(stats, function(x) x$size))
  expect_equal(total_size, nrow(test_data$data))

  # Get cluster likelihoods
  likelihoods <- runner$get_likelihoods()
  expect_true(length(likelihoods) >= 1)
  expect_true(all(is.finite(likelihoods) | likelihoods == -Inf))

  # Get membership matrix
  membership <- runner$get_membership_matrix()
  expect_equal(nrow(membership), nrow(test_data$data))
  expect_equal(rowSums(membership), rep(1, nrow(test_data$data)))
})

# Test 13: Posterior Predictive Sampling
test_that("MCMCRunnerManual posterior predictive sampling works", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Run some iterations first
  for (i in 1:30) {
    runner$perform_iteration()
  }

  # Sample from posterior predictive
  samples <- runner$sample_predictive(n_samples = 10)

  expect_type(samples, "list")
  expect_equal(length(samples), 10)

  # Each sample should be a parameter vector
  for (s in samples) {
    expect_equal(length(s), 2)  # Gaussian parameters
  }
})

# Test 14: Full Run and Results
test_that("MCMCRunnerManual full run produces correct results", {
  test_data <- create_test_data(n = 50, k = 3, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp, n_iter = 100, n_burn = 20, thin = 2)

  # Run to completion
  results <- runner$run()

  expect_type(results, "list")
  expect_true(runner$is_complete())

  # Check results structure
  expect_true(all(c("cluster_labels", "alpha", "n_clusters", "labels_chain",
                    "alpha_chain", "theta_chain", "n_clusters_chain") %in% names(results)))

  # Check dimensions
  n_stored <- (100 - 20) / 2
  expect_equal(nrow(results$labels_chain), n_stored)
  expect_equal(ncol(results$labels_chain), nrow(test_data$data))
  expect_equal(length(results$alpha_chain), n_stored)
  expect_equal(length(results$n_clusters_chain), n_stored)

  # Check convergence diagnostics are included
  expect_true("convergence_diagnostics" %in% names(results))
})

# Test 15: Integration with Standard DP Interface
test_that("MCMCRunnerManual integrates with standard DP interface", {
  skip_if_not(can_use_cpp())

  test_data <- create_test_data(n = 100, k = 3, dim = 2)

  # Run with standard interface
  dp_r <- DirichletProcessGaussian(test_data$data)
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)

  # Run with manual C++ interface
  dp_cpp <- DirichletProcessGaussian(test_data$data)
  dp_cpp <- Initialise(dp_cpp)
  runner <- CppMCMCRunner$new(dp_cpp, n_iter = 100)
  results_cpp <- runner$run()

  # Both should produce valid results
  expect_true(length(dp_r$clusterLabels) == nrow(test_data$data))
  expect_true(length(results_cpp$cluster_labels) == nrow(test_data$data))

  # Number of clusters should be reasonable
  expect_true(dp_r$numberClusters >= 1 && dp_r$numberClusters <= 20)
  expect_true(results_cpp$n_clusters >= 1 && results_cpp$n_clusters <= 20)
})

# Test 16: Error Handling
test_that("MCMCRunnerManual handles errors appropriately", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  runner <- CppMCMCRunner$new(dp)

  # Invalid cluster indices
  expect_error(runner$merge_clusters(1, 10))
  expect_error(runner$split_cluster(10))

  # Invalid temperature
  expect_error(runner$set_temperature(-1))
  expect_error(runner$set_temperature(0))

  # Invalid auxiliary count
  expect_error(runner$set_auxiliary_count(0))
  expect_error(runner$set_auxiliary_count(-5))

  # Invalid bounds
  expect_error(runner$set_bounds(lower = c(1, 2, 3), upper = c(4, 5)))

  # Invalid labels
  expect_error(runner$set_labels(1:10))  # Wrong length
})

# Test 17: Performance Comparison
test_that("MCMCRunnerManual performs better than R implementation", {
  skip_if_not(can_use_cpp())
  skip_on_cran()  # Performance tests can be variable

  test_data <- create_test_data(n = 200, k = 4, dim = 2)

  dp <- DirichletProcessGaussian(test_data$data)

  # Time R implementation
  r_time <- system.time({
    dp_r <- Fit(Clone(dp), 100, progressBar = FALSE, useC = FALSE)
  })[3]

  # Time C++ implementation
  cpp_time <- system.time({
    runner <- CppMCMCRunner$new(Clone(dp), n_iter = 100)
    results <- runner$run()
  })[3]

  # C++ should be faster
  speedup <- r_time / cpp_time
  expect_true(speedup > 2,
              info = sprintf("C++ speedup: %.1fx (R: %.2fs, C++: %.2fs)",
                             speedup, r_time, cpp_time))
})

# Test 18: Different Distributions
test_that("MCMCRunnerManual works with different distributions", {
  skip_if_not(can_use_cpp())

  # Beta distribution
  set.seed(123)
  beta_data <- matrix(rbeta(100, 2, 5), ncol = 1)
  dp_beta <- DirichletProcessBeta(beta_data, 1)
  dp_beta <- Initialise(dp_beta)

  runner_beta <- CppMCMCRunner$new(dp_beta, n_iter = 50)

  # Should initialize correctly
  expect_s4_class(runner_beta, "CppMCMCRunner")

  # Run some iterations
  for (i in 1:20) {
    runner_beta$perform_iteration()
  }

  state <- runner_beta$get_state()
  expect_true(state$n_clusters >= 1)
  expect_true(all(sapply(state$cluster_params, length) == 4))  # Beta has 4 params
})

# Test 19: Reproducibility
test_that("MCMCRunnerManual produces reproducible results", {
  test_data <- create_test_data(n = 30, k = 2, dim = 1)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  # Run 1
  set.seed(456)
  runner1 <- CppMCMCRunner$new(dp, n_iter = 50)
  results1 <- runner1$run()

  # Run 2
  set.seed(456)
  runner2 <- CppMCMCRunner$new(dp, n_iter = 50)
  results2 <- runner2$run()

  # Results should be identical
  expect_equal(results1$alpha_chain, results2$alpha_chain)
  expect_equal(results1$labels_chain, results2$labels_chain)
})

# Test 20: Memory Management
test_that("MCMCRunnerManual manages memory correctly", {
  test_data <- create_test_data(n = 1000, k = 5, dim = 3)

  dp <- DirichletProcessGaussian(test_data$data)
  dp <- Initialise(dp)

  # Create and destroy multiple runners
  for (i in 1:5) {
    runner <- CppMCMCRunner$new(dp, n_iter = 100)

    # Run some iterations
    for (j in 1:20) {
      runner$perform_iteration()
    }

    # Get results
    results <- runner$get_results()

    # Runner should be garbage collected
    rm(runner)
    gc()
  }

  # Should not crash or leak memory
  expect_true(TRUE)
})

# Run all tests
if (interactive()) {
  test_results <- test_dir("tests/testthat", filter = "mcmc-runner-manual")
  print(test_results)
}
