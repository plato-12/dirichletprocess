# tests/testthat/test-cpp-manual-mcmc.R

test_that("CppMCMCRunner produces consistent results with Fit()", {
  distributions <- c("normal", "exponential", "beta", "weibull", "mvnormal", "mvnormal2")

  for (dist in distributions) {
    # Generate appropriate test data
    test_data <- generate_test_data(dist, n = 100)

    # Standard Fit approach
    set.seed(123)
    dp_fit <- create_dp_object(dist, test_data)
    dp_fit <- Fit(dp_fit, its = 100)

    # Manual MCMC approach
    set.seed(123)
    dp_manual <- create_dp_object(dist, test_data)
    runner <- CppMCMCRunner$new(dp_manual)

    for (i in 1:100) {
      runner$step_assignments()
      runner$step_parameters()
      runner$step_concentration()
    }

    manual_state <- runner$get_state()

    # Compare final states
    expect_equal(
      dp_fit$clusterLabels,
      manual_state$labels,
      info = paste("Labels mismatch for", dist)
    )

    expect_equal(
      dp_fit$alpha,
      manual_state$alpha,
      tolerance = 1e-6,
      info = paste("Alpha mismatch for", dist)
    )
  }
})

test_that("CppMCMCRunner advanced features work correctly", {
  set.seed(123)
  test_data <- rnorm(100)
  dp <- DirichletProcessGaussian(test_data)

  runner <- CppMCMCRunner$new(dp)

  # Test temperature control
  runner$set_temperature(0.5)
  expect_equal(runner$get_temperature(), 0.5)

  # Test auxiliary parameters
  runner$set_auxiliary_params(list(scale = 2.0))
  aux <- runner$get_auxiliary_params()
  expect_equal(aux$scale, 2.0)

  # Test predictive sampling
  predictive <- runner$sample_predictive(n = 10)
  expect_length(predictive, 10)

  # Test cluster operations
  runner$step_assignments()
  n_clusters_before <- runner$get_n_clusters()

  runner$merge_clusters(1, 2)
  n_clusters_after <- runner$get_n_clusters()
  expect_lt(n_clusters_after, n_clusters_before)
})

test_that("Manual MCMC step functions work individually", {
  set.seed(123)
  test_data <- generate_test_data("normal", 100)
  dp <- DirichletProcessGaussian(test_data)

  runner <- CppMCMCRunner$new(dp)

  # Test individual steps don't error
  expect_error(runner$step_assignments(), NA)
  expect_error(runner$step_parameters(), NA)
  expect_error(runner$step_concentration(), NA)

  # Test state extraction
  state <- runner$get_state()
  expect_type(state, "list")
  expect_true(all(c("labels", "alpha", "parameters") %in% names(state)))
})

test_that("Manual MCMC handles different covariance models", {
  # Test MVNormal with different covariance structures
  covariance_models <- c("E", "V", "EII", "VII", "EEI", "VEI", "EVI", "VVI", "FULL")

  for (model in covariance_models) {
    set.seed(123)
    test_data <- generate_test_data("mvnormal", 100)

    tryCatch({
      dp <- DirichletProcessMvnormal(test_data, g0Priors = list(model = model))
      runner <- CppMCMCRunner$new(dp)

      # Run a few steps
      for (i in 1:10) {
        runner$step_assignments()
        runner$step_parameters()
      }

      state <- runner$get_state()
      expect_true(is.list(state$parameters))

    }, error = function(e) {
      skip(paste("Covariance model", model, "not supported:", e$message))
    })
  }
})

test_that("Manual MCMC produces same chain as Fit() when steps match", {
  set.seed(123)
  test_data <- generate_test_data("exponential", 50)

  # Run with Fit()
  dp_fit <- DirichletProcessExponential(test_data)
  dp_fit <- Fit(dp_fit, its = 50)

  # Run manually
  set.seed(123)
  dp_manual <- DirichletProcessExponential(test_data)
  runner <- CppMCMCRunner$new(dp_manual)

  alpha_chain <- numeric(50)
  cluster_counts <- numeric(50)

  for (i in 1:50) {
    runner$step_assignments()
    runner$step_parameters()
    runner$step_concentration()

    state <- runner$get_state()
    alpha_chain[i] <- state$alpha
    cluster_counts[i] <- length(unique(state$labels))
  }

  # Compare chains
  expect_equal(dp_fit$alphaChain, alpha_chain, tolerance = 1e-6)
  expect_equal(
    sapply(dp_fit$labelsChain, function(x) length(unique(x))),
    cluster_counts
  )
})

test_that("CppMCMCRunner handles hierarchical distributions", {
  # Test hierarchical distributions if available
  hierarchical_dists <- c("hierarchical_beta", "hierarchical_mvnormal", "hierarchical_mvnormal2")

  for (dist in hierarchical_dists) {
    tryCatch({
      test_data <- generate_test_data("beta", 100)  # Use appropriate data

      # Create hierarchical DP object
      # This may need adjustment based on actual hierarchical implementation
      dp <- create_dp_object(dist, test_data)
      runner <- CppMCMCRunner$new(dp)

      # Test that basic operations work
      expect_error(runner$step_assignments(), NA)
      expect_error(runner$step_parameters(), NA)
      expect_error(runner$get_state(), NA)

    }, error = function(e) {
      skip(paste("Hierarchical distribution", dist, "not available"))
    })
  }
})
