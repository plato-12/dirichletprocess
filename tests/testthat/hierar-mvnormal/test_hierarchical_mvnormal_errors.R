context("Hierarchical MVNormal Error Handling Tests")

test_that("Hierarchical MVNormal handles invalid inputs gracefully", {
  # Test with invalid prior parameters
  expect_error(
    hierarchical_mvnormal_create_mixing(
      n_groups = 2,
      prior_params = list(
        mu0 = c(0, 0),
        kappa0 = -1,  # Invalid: negative kappa
        Lambda = diag(2),
        nu = 4
      ),
      alpha_prior = c(2, 2),
      gamma_prior = c(2, 2),
      n_sticks = 10
    )
  )

  # Test with mismatched dimensions
  expect_error(
    hierarchical_mvnormal_create_mixing(
      n_groups = 2,
      prior_params = list(
        mu0 = c(0, 0, 0),  # 3D
        kappa0 = 1,
        Lambda = diag(2),  # 2D - mismatch
        nu = 4
      ),
      alpha_prior = c(2, 2),
      gamma_prior = c(2, 2),
      n_sticks = 10
    )
  )

  # Test with invalid number of groups
  expect_error(
    hierarchical_mvnormal_create_mixing(
      n_groups = 0,  # Invalid
      prior_params = list(mu0 = 0, kappa0 = 1, Lambda = matrix(1), nu = 2),
      alpha_prior = c(2, 2),
      gamma_prior = c(2, 2),
      n_sticks = 10
    )
  )
})

test_that("Hierarchical MVNormal MCMC handles edge cases", {
  # Empty data list
  expect_error(
    hierarchical_mvnormal_run(
      data_list = list(),
      hdp_params = list(
        n_sticks = 10,
        prior_params = list(mu0 = 0, kappa0 = 1, Lambda = matrix(1), nu = 2),
        alpha_prior = c(2, 2),
        gamma_prior = c(2, 2)
      ),
      mcmc_params = list(n_iter = 10, n_burn = 0, thin = 1,
                         update_prior = TRUE, show_progress = FALSE)
    )
  )

  # Mismatched data dimensions
  data_list_bad <- list(
    matrix(rnorm(20), ncol = 2),
    matrix(rnorm(30), ncol = 3)  # Different dimension
  )

  expect_error(
    hierarchical_mvnormal_run(
      data_list_bad,
      hdp_params = list(
        n_sticks = 10,
        prior_params = list(mu0 = c(0,0), kappa0 = 1, Lambda = diag(2), nu = 4),
        alpha_prior = c(2, 2),
        gamma_prior = c(2, 2)
      ),
      mcmc_params = list(n_iter = 10, n_burn = 0, thin = 1,
                         update_prior = TRUE, show_progress = FALSE)
    )
  )
})

test_that("Hierarchical MVNormal handles numerical edge cases", {
  # Very small variance data
  set.seed(222)
  data_list <- list(
    matrix(rep(c(1, 2), 20), ncol = 2, byrow = TRUE) +
      matrix(rnorm(40, 0, 0.001), ncol = 2),
    matrix(rep(c(-1, -2), 20), ncol = 2, byrow = TRUE) +
      matrix(rnorm(40, 0, 0.001), ncol = 2)
  )

  hdp_params <- list(
    n_sticks = 5,
    prior_params = list(
      mu0 = c(0, 0),
      kappa0 = 0.01,  # Small kappa
      Lambda = diag(2) * 100,  # Large precision
      nu = 3
    ),
    alpha_prior = c(0.1, 0.1),  # Small alphas
    gamma_prior = c(0.1, 0.1)   # Small gamma prior
  )

  mcmc_params <- list(
    n_iter = 20,
    n_burn = 0,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  # Should handle without numerical issues
  expect_error(
    result <- hierarchical_mvnormal_run(data_list, hdp_params, mcmc_params),
    NA
  )

  # Check results are finite
  final_state <- result$final_state
  expect_true(all(is.finite(final_state$stick_weights)))
  expect_true(all(is.finite(final_state$alphas)))
  expect_true(is.finite(final_state$gamma))
})
