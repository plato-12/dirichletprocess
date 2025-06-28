context("Hierarchical MVNormal Integration Tests")

test_that("Hierarchical MVNormal integrates with R interface", {
  skip_if_not(using_cpp())
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Create hierarchical data
  set.seed(789)
  dataList <- list(
    mvtnorm::rmvnorm(30, c(-1, -1), diag(2)),
    mvtnorm::rmvnorm(30, c(1, 1), diag(2)),
    mvtnorm::rmvnorm(30, c(0, 2), diag(2))
  )

  # Use the R interface (if it exists)
  # This tests the full integration pathway
  g0Priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2) * 2,
    kappa0 = 0.5,
    nu = 4
  )

  # Direct C++ call
  hdp_params <- list(
    n_sticks = 15,
    prior_params = g0Priors,
    alpha_prior = c(2, 4),
    gamma_prior = c(2, 2)
  )

  mcmc_params <- list(
    n_iter = 50,
    n_burn = 10,
    thin = 2,
    update_prior = TRUE,
    show_progress = FALSE
  )

  result <- hierarchical_mvnormal_run(dataList, hdp_params, mcmc_params)

  # Verify output can be used for posterior inference
  expect_true(length(result$samples) == 20)  # (50-10)/2

  # Extract posterior samples of cluster parameters
  final_params <- result$samples[[20]]$cluster_params
  expect_equal(length(final_params), 3)  # 3 groups

  # Each group should have reasonable parameter dimensions
  for (g in 1:3) {
    if (length(final_params[[g]]) > 0) {
      # Each parameter should be a 2D vector
      expect_equal(length(final_params[[g]][[1]]), 2)
    }
  }
})

test_that("Hierarchical MVNormal posterior sampling works", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Generate data
  dataList <- list(
    mvtnorm::rmvnorm(20, c(0, 0), diag(2)),
    mvtnorm::rmvnorm(20, c(0, 0), diag(2))
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
    gamma_prior = c(2, 2)
  )

  mcmc_params <- list(
    n_iter = 100,
    n_burn = 50,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  # Fit model
  result <- hierarchical_mvnormal_run(dataList, hdp_params, mcmc_params)

  # Test posterior sampling functionality
  n_post_samples <- 100
  post_samples <- hierarchical_mvnormal_posterior_sample(
    result$final_state,
    dataList,
    n_post_samples
  )

  expect_type(post_samples, "list")
  expect_equal(length(post_samples), n_post_samples)

  # Each sample should contain parameters for each group
  for (i in 1:10) {  # Check first 10 samples
    expect_equal(length(post_samples[[i]]), 2)  # 2 groups
  }
})

test_that("Hierarchical MVNormal preserves data integrity", {
  # Ensure data is not modified during MCMC
  set.seed(111)

  original_data <- list(
    matrix(rnorm(40), ncol = 2),
    matrix(rnorm(40), ncol = 2)
  )

  # Make a deep copy
  data_list <- list(
    original_data[[1]][,],
    original_data[[2]][,]
  )

  hdp_params <- list(
    n_sticks = 5,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 1,
      Lambda = diag(2),
      nu = 4
    ),
    alpha_prior = c(2, 2),
    gamma_prior = c(2, 2)
  )

  mcmc_params <- list(
    n_iter = 20,
    n_burn = 0,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params)

  # Check data integrity
  expect_equal(result$data[[1]], original_data[[1]])
  expect_equal(result$data[[2]], original_data[[2]])
})
