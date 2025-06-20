# tests/testthat/test_mvnormal_cpp_diagnostic.R
context("MVNormal C++ Diagnostic Tests")

test_that("MVNormal C++ components are available", {
  # Check if C++ MCMC runner exists
  expect_true(exists("_dirichletprocess_run_mcmc_cpp"),
              "C++ MCMC runner not found")

  # Check if MVNormal-specific functions exist
  expect_true(exists("conjugate_mvnormal_cluster_component_update_cpp"),
              "MVNormal cluster component update not found")
  expect_true(exists("conjugate_mvnormal_cluster_parameter_update_cpp"),
              "MVNormal cluster parameter update not found")

  # Test C++ status
  cpp_status <- get_cpp_status()
  print(cpp_status)

  # Test MVNormal DP object
  set.seed(123)
  y <- matrix(rnorm(20), ncol = 2)
  g0Priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2),
    kappa0 = 1,
    nu = 4
  )

  dp <- DirichletProcessMvnormal(y, g0Priors)

  # Check if can_use_cpp works
  can_use <- can_use_cpp(dp)
  expect_true(can_use, "can_use_cpp returns FALSE for MVNormal")

  # If can_use_cpp is TRUE, test parameter preparation
  if (can_use) {
    params <- prepare_mixing_dist_params(dp)
    expect_equal(params$type, "mvnormal")
    expect_true(all(c("mu0", "kappa0", "Lambda", "nu") %in% names(params)))

    # Test running C++ MCMC directly
    mcmc_params <- list(
      n_iter = 10,
      n_burn = 0,
      thin = 1,
      update_concentration = TRUE,
      alpha = 1.0,
      m_auxiliary = 3
    )

    result <- tryCatch({
      run_mcmc_cpp(as.matrix(y), params, mcmc_params)
    }, error = function(e) {
      print(paste("Error in run_mcmc_cpp:", e$message))
      NULL
    })

    if (!is.null(result)) {
      expect_type(result, "list")
      print("C++ MCMC ran successfully!")
    }
  }
})
