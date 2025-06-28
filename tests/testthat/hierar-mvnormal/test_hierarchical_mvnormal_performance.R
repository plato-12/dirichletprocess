context("Hierarchical MVNormal C++ Performance Tests")

test_that("Hierarchical MVNormal scales with number of groups", {
  skip_on_cran()  # Skip on CRAN due to time constraints
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Test scaling with different numbers of groups
  n_groups_vec <- c(2, 5, 10)
  times <- numeric(length(n_groups_vec))

  for (i in seq_along(n_groups_vec)) {
    n_groups <- n_groups_vec[i]

    # Generate data
    data_list <- list()
    for (g in 1:n_groups) {
      data_list[[g]] <- matrix(rnorm(100), ncol = 2)
    }

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
      n_iter = 50,
      n_burn = 0,
      thin = 1,
      update_prior = TRUE,
      show_progress = FALSE
    )

    # Time the execution
    start_time <- Sys.time()
    result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params)
    end_time <- Sys.time()

    times[i] <- as.numeric(end_time - start_time, units = "secs")
  }

  # Check that scaling is reasonable (not quadratic)
  # Time should increase roughly linearly with number of groups
  time_ratio <- times[3] / times[1]
  group_ratio <- n_groups_vec[3] / n_groups_vec[1]

  # Allow for some overhead, but should be less than quadratic
  expect_true(time_ratio < group_ratio^1.5)
})

test_that("Hierarchical MVNormal handles large datasets efficiently", {
  skip_on_cran()
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Test with larger dataset
  n_obs <- 500
  n_groups <- 3

  data_list <- list()
  for (g in 1:n_groups) {
    data_list[[g]] <- matrix(rnorm(n_obs * 2), ncol = 2)
  }

  hdp_params <- list(
    n_sticks = 20,
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

  # Should complete in reasonable time
  start_time <- Sys.time()
  result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params)
  end_time <- Sys.time()

  time_taken <- as.numeric(end_time - start_time, units = "secs")

  # Should complete within 30 seconds even for large data
  expect_true(time_taken < 30)

  # Check that result is valid
  expect_equal(length(result$samples), 20)
})
