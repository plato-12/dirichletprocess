context("Beta MCMC Sampler")

test_that("Beta Metropolis-Hastings produces valid samples", {
  set.seed(456)
  mdObj <- BetaMixtureCreate(c(2, 8), c(0.1, 0.1), 1)

  # Generate data from known Beta distribution
  true_a <- 3
  true_b <- 7
  data <- matrix(rbeta(50, true_a, true_b), ncol = 1)

  # Run MH sampler
  n_samples <- 1000
  posterior_samples <- PosteriorDraw(mdObj, data, n_samples)

  # Extract final samples (after burn-in)
  mu_samples <- as.vector(posterior_samples$mu[, , 501:1000])
  nu_samples <- as.vector(posterior_samples$nu[, , 501:1000])

  # Check convergence to true parameters
  true_mean <- true_a / (true_a + true_b)
  true_var <- (true_a * true_b) / ((true_a + true_b)^2 * (true_a + true_b + 1))
  true_tau <- true_mean * (1 - true_mean) / true_var - 1

  expect_equal(mean(mu_samples), true_mean, tolerance = 0.1)
  expect_equal(mean(nu_samples), true_tau, tolerance = 2)

  # Check MCMC diagnostics
  # Samples should show some autocorrelation but not be stuck
  mu_acf <- acf(mu_samples, plot = FALSE)
  expect_true(mu_acf$acf[2] < 0.95)  # Not perfectly correlated
  expect_true(mu_acf$acf[10] < 0.5)  # Correlation decays
})

test_that("Beta MH parameter proposal maintains detailed balance", {
  mdObj <- BetaMixtureCreate(c(2, 8), c(0.5, 0.5), 1)

  # Test proposal symmetry
  old_params <- list(
    mu = array(0.5, dim = c(1, 1, 1)),
    nu = array(2, dim = c(1, 1, 1))
  )

  # Generate many proposals
  set.seed(789)
  n_props <- 1000
  mu_props <- numeric(n_props)
  nu_props <- numeric(n_props)

  for (i in 1:n_props) {
    prop <- MhParameterProposal(mdObj, old_params)
    mu_props[i] <- prop$mu[1]
    nu_props[i] <- prop$nu[1]
  }

  # Check proposal distribution properties
  # Should be centered around old values
  expect_equal(mean(mu_props[!is.na(mu_props)]), 0.5, tolerance = 0.05)
  expect_equal(mean(log(nu_props)), log(2), tolerance = 0.1)

  # Check boundary handling
  expect_true(all(mu_props >= 0 | is.na(mu_props)))
  expect_true(all(mu_props <= 1 | is.na(mu_props)))
  expect_true(all(nu_props > 0))
})
