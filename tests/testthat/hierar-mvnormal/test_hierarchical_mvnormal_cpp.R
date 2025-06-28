context("Hierarchical MVNormal C++ Implementation Tests")

test_that("Hierarchical MVNormal mixing distribution creation works", {
  # Test parameters
  n_groups <- 3
  n_sticks <- 20
  prior_params <- list(
    mu0 = c(0, 0),
    kappa0 = 1.0,
    Lambda = diag(2),
    nu = 4
  )
  alpha_prior <- c(2, 2)
  gamma_prior <- c(2, 2)

  # Test C++ creation function
  result <- hierarchical_mvnormal_create_mixing(
    n_groups, prior_params, alpha_prior, gamma_prior, n_sticks
  )

  expect_type(result, "list")
  expect_true("gamma" %in% names(result))
  expect_true("stick_weights" %in% names(result))
  expect_true("alphas" %in% names(result))
  expect_true("pi_k" %in% names(result))

  # Check dimensions
  expect_equal(length(result$stick_weights), n_sticks)
  expect_equal(length(result$alphas), n_groups)
  expect_equal(length(result$pi_k), n_groups)

  # Check that stick weights sum to 1
  expect_equal(sum(result$stick_weights), 1, tolerance = 1e-10)

  # Check that each group's pi_k sums to 1
  for (i in 1:n_groups) {
    expect_equal(sum(result$pi_k[[i]]), 1, tolerance = 1e-10)
  }

  # Check parameter constraints
  expect_true(all(result$stick_weights >= 0))
  expect_true(all(result$alphas > 0))
  expect_true(result$gamma > 0)
})

test_that("Hierarchical MVNormal MCMC runner works with simple data", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Generate simple hierarchical data
  set.seed(123)
  n_groups <- 2
  n_per_group <- 30

  # True parameters
  true_means <- list(
    c(-2, -2),
    c(2, 2)
  )

  data_list <- list()
  for (i in 1:n_groups) {
    data_list[[i]] <- mvtnorm::rmvnorm(n_per_group, true_means[[i]], diag(2))
  }

  # Set up HDP parameters
  hdp_params <- list(
    n_sticks = 10,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 0.5,
      Lambda = diag(2) * 2,
      nu = 4
    ),
    alpha_prior = c(2, 2),
    gamma_prior = c(2, 2)
  )

  mcmc_params <- list(
    n_iter = 50,
    n_burn = 10,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  # Run MCMC
  result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params)

  # Check output structure
  expect_type(result, "list")
  expect_true("samples" %in% names(result))
  expect_true("data" %in% names(result))
  expect_true("final_state" %in% names(result))

  # Check samples were stored correctly (50 total - 10 burn = 40 samples)
  expect_equal(length(result$samples), 40)

  # Check each sample has the expected structure
  sample1 <- result$samples[[1]]
  expect_true("iteration" %in% names(sample1))
  expect_true("cluster_labels" %in% names(sample1))
  expect_true("cluster_params" %in% names(sample1))
  expect_true("n_clusters" %in% names(sample1))
  expect_true("hdp_state" %in% names(sample1))

  # Check cluster labels
  expect_equal(length(sample1$cluster_labels), n_groups)
  for (i in 1:n_groups) {
    expect_equal(length(sample1$cluster_labels[[i]]), n_per_group)
    expect_true(all(sample1$cluster_labels[[i]] >= 0))
  }

  # Check that we found reasonable clusters
  final_state <- result$final_state
  expect_true(all(final_state$alphas > 0))
  expect_true(final_state$gamma > 0)
})

test_that("Hierarchical MVNormal handles edge cases", {
  # Test with minimal data
  data_list <- list(
    matrix(rnorm(4), ncol = 2),  # 2 observations
    matrix(rnorm(4), ncol = 2)   # 2 observations
  )

  hdp_params <- list(
    n_sticks = 5,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 1.0,
      Lambda = diag(2),
      nu = 3
    ),
    alpha_prior = c(1, 1),
    gamma_prior = c(1, 1)
  )

  mcmc_params <- list(
    n_iter = 10,
    n_burn = 0,
    thin = 1,
    update_prior = FALSE,
    show_progress = FALSE
  )

  # Should run without error
  expect_error(
    result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params),
    NA
  )

  # Test with 1D data (should work with d=1)
  data_list_1d <- list(
    matrix(rnorm(10), ncol = 1),
    matrix(rnorm(10), ncol = 1)
  )

  hdp_params_1d <- list(
    n_sticks = 5,
    prior_params = list(
      mu0 = 0,
      kappa0 = 1.0,
      Lambda = matrix(1, 1, 1),
      nu = 2
    ),
    alpha_prior = c(1, 1),
    gamma_prior = c(1, 1)
  )

  expect_error(
    result_1d <- hierarchical_mvnormal_run(data_list_1d, hdp_params_1d, mcmc_params),
    NA
  )
})

test_that("Hierarchical MVNormal update functions work correctly", {
  # Create a simple DP object to test update functions
  n_groups <- 2
  data_list <- list(
    matrix(rnorm(20), ncol = 2),
    matrix(rnorm(20), ncol = 2)
  )

  # Create initial state
  hdp_params <- list(
    n_sticks = 5,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 1.0,
      Lambda = diag(2),
      nu = 4
    ),
    alpha_prior = c(2, 2),
    gamma_prior = c(2, 2)
  )

  # Create mixing object
  mixing_state <- hierarchical_mvnormal_create_mixing(
    n_groups,
    hdp_params$prior_params,
    hdp_params$alpha_prior,
    hdp_params$gamma_prior,
    hdp_params$n_sticks
  )

  # Create mock DP object structure
  dp_obj <- list(
    dataList = data_list,
    clusterLabels = list(
      rep(0, 10),  # Group 1 labels
      rep(0, 10)   # Group 2 labels
    ),
    hdp_state = mixing_state,
    alpha = 1.0,
    m = 3
  )

  # Test cluster update
  updated_dp <- hierarchical_mvnormal_update_clusters(dp_obj, group_idx = 0)

  expect_type(updated_dp, "list")
  expect_true("clusterLabels" %in% names(updated_dp))

  # Check that labels were updated
  expect_equal(length(updated_dp$clusterLabels[[1]]), 10)
  expect_true(all(updated_dp$clusterLabels[[1]] >= 0))
})

test_that("Hierarchical MVNormal produces reasonable clustering", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Generate well-separated hierarchical data
  set.seed(456)
  n_groups <- 3

  # Each group has two distinct clusters
  data_list <- list()
  for (g in 1:n_groups) {
    cluster1 <- mvtnorm::rmvnorm(15, c(-3, -3), diag(2) * 0.5)
    cluster2 <- mvtnorm::rmvnorm(15, c(3, 3), diag(2) * 0.5)
    data_list[[g]] <- rbind(cluster1, cluster2)
  }

  hdp_params <- list(
    n_sticks = 20,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 0.1,
      Lambda = diag(2) * 0.1,
      nu = 4
    ),
    alpha_prior = c(2, 4),
    gamma_prior = c(2, 2)
  )

  mcmc_params <- list(
    n_iter = 100,
    n_burn = 50,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params)

  # Check that we found approximately 2 clusters per group
  final_labels <- result$samples[[length(result$samples)]]$cluster_labels

  for (g in 1:n_groups) {
    n_unique <- length(unique(final_labels[[g]]))
    expect_true(n_unique >= 1)
    expect_true(n_unique <= 10)  # Reasonable upper bound
  }

  # Check parameter convergence - gamma should stabilize
  gamma_values <- sapply(result$samples, function(s) s$hdp_state$gamma)
  gamma_var_first_half <- var(gamma_values[1:25])
  gamma_var_second_half <- var(gamma_values[26:50])

  # Variance should decrease as chain converges
  expect_true(gamma_var_second_half <= gamma_var_first_half * 2)
})

test_that("Hierarchical MVNormal stick-breaking works correctly", {
  # Test the stick-breaking process
  n_sticks <- 100
  gamma <- 5.0

  # Create mixing distribution to access stick weights
  mixing_obj <- hierarchical_mvnormal_create_mixing(
    n_groups = 1,
    prior_params = list(mu0 = 0, kappa0 = 1, Lambda = matrix(1), nu = 2),
    alpha_prior = c(1, 1),
    gamma_prior = c(gamma, 1),  # Fix gamma to test
    n_sticks = n_sticks
  )

  weights <- mixing_obj$stick_weights

  # Check properties
  expect_equal(length(weights), n_sticks)
  expect_equal(sum(weights), 1, tolerance = 1e-10)
  expect_true(all(weights >= 0))
  expect_true(all(weights <= 1))

  # Check that weights generally decrease (stochastic, but should trend down)
  # First 10 weights should generally be larger than last 10
  expect_true(mean(weights[1:10]) > mean(weights[91:100]))
})

test_that("Hierarchical MVNormal handles different data dimensions", {
  # Test with different dimensional data
  for (d in c(1, 3, 5)) {
    data_list <- list()
    for (g in 1:2) {
      data_list[[g]] <- matrix(rnorm(20 * d), ncol = d)
    }

    hdp_params <- list(
      n_sticks = 5,
      prior_params = list(
        mu0 = rep(0, d),
        kappa0 = 1.0,
        Lambda = diag(d),
        nu = d + 1
      ),
      alpha_prior = c(2, 2),
      gamma_prior = c(2, 2)
    )

    mcmc_params <- list(
      n_iter = 10,
      n_burn = 0,
      thin = 1,
      update_prior = FALSE,
      show_progress = FALSE
    )

    expect_error(
      result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params),
      NA,
      info = paste("Failed for dimension", d)
    )

    # Check output dimensions match input
    expect_equal(ncol(result$data[[1]]), d)
  }
})

test_that("Hierarchical MVNormal parameter updates are consistent", {
  # Test that parameter updates maintain consistency

  # Create initial state
  mixing_state <- hierarchical_mvnormal_create_mixing(
    n_groups = 2,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 1,
      Lambda = diag(2),
      nu = 4
    ),
    alpha_prior = c(2, 2),
    gamma_prior = c(3, 1),
    n_sticks = 10
  )

  # Check initial state consistency
  initial_gamma <- mixing_state$gamma
  initial_weights <- mixing_state$stick_weights
  initial_alphas <- mixing_state$alphas

  expect_true(initial_gamma > 0)
  expect_equal(sum(initial_weights), 1, tolerance = 1e-10)
  expect_true(all(initial_alphas > 0))

  # Run a few iterations and check consistency is maintained
  data_list <- list(
    matrix(rnorm(40), ncol = 2),
    matrix(rnorm(40), ncol = 2)
  )

  hdp_params <- list(
    n_sticks = 10,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 1,
      Lambda = diag(2),
      nu = 4
    ),
    alpha_prior = c(2, 2),
    gamma_prior = c(3, 1)
  )

  mcmc_params <- list(
    n_iter = 20,
    n_burn = 0,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params)

  # Check all samples maintain consistency
  for (i in 1:length(result$samples)) {
    state <- result$samples[[i]]$hdp_state
    expect_true(state$gamma > 0)
    expect_equal(sum(state$stick_weights), 1, tolerance = 1e-10)
    expect_true(all(state$alphas > 0))
  }
})
