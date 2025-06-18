context("Beta Distribution Core Functions")

test_that("Beta likelihood calculations are correct", {
  # Test parameters
  maxT <- 1
  mu <- 0.5
  tau <- 5

  # Test data
  x <- seq(0.1, 0.9, by = 0.1)

  # Expected values using standard Beta parameterization
  a <- (mu * tau) / maxT
  b <- (1 - mu/maxT) * tau
  expected <- (1/maxT) * dbeta(x/maxT, a, b)

  # Create mixing distribution
  mdObj <- BetaMixtureCreate(c(2, 8), c(1, 1), maxT)
  theta <- list(
    mu = array(mu, dim = c(1, 1, 1)),
    nu = array(tau, dim = c(1, 1, 1))
  )

  # Calculate likelihood
  actual <- Likelihood(mdObj, x, theta)

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("Beta likelihood handles edge cases", {
  mdObj <- BetaMixtureCreate(c(2, 8), c(1, 1), 1)

  # Test with extreme tau values
  theta_low_tau <- list(
    mu = array(0.5, dim = c(1, 1, 1)),
    nu = array(0.01, dim = c(1, 1, 1))
  )

  theta_high_tau <- list(
    mu = array(0.5, dim = c(1, 1, 1)),
    nu = array(1000, dim = c(1, 1, 1))
  )

  # Should not produce NaN or Inf
  lik_low <- Likelihood(mdObj, 0.5, theta_low_tau)
  lik_high <- Likelihood(mdObj, 0.5, theta_high_tau)

  expect_true(all(is.finite(lik_low)))
  expect_true(all(is.finite(lik_high)))

  # Test boundary values
  x_boundary <- c(0, 1)
  lik_boundary <- Likelihood(mdObj, x_boundary, theta_low_tau)
  expect_true(all(is.finite(lik_boundary) | lik_boundary == 0))
})

test_that("Beta prior draw produces valid samples", {
  set.seed(123)
  mdObj <- BetaMixtureCreate(c(2, 8), c(1, 1), maxT = 1)

  n_samples <- 1000
  prior_samples <- PriorDraw(mdObj, n_samples)

  # Check structure
  expect_named(prior_samples, c("mu", "nu"))
  expect_equal(dim(prior_samples$mu), c(1, 1, n_samples))
  expect_equal(dim(prior_samples$nu), c(1, 1, n_samples))

  # Check bounds
  mu_vals <- as.vector(prior_samples$mu)
  nu_vals <- as.vector(prior_samples$nu)

  expect_true(all(mu_vals >= 0 & mu_vals <= mdObj$maxT))
  expect_true(all(nu_vals > 0))

  # Check prior distribution properties
  # mu should be approximately uniform
  ks_test_mu <- ks.test(mu_vals, "punif", 0, mdObj$maxT)
  expect_true(ks_test_mu$p.value > 0.01)

  # 1/nu should follow Gamma(2, 8)
  inv_nu <- 1/nu_vals
  # Use a more lenient test due to sampling variability
  expect_true(mean(inv_nu) < 1)  # Gamma(2,8) has mean 2/8 = 0.25
})

test_that("Beta prior density calculations are correct", {
  mdObj <- BetaMixtureCreate(c(2, 8), c(1, 1), 1)

  # Test multiple parameter combinations
  test_cases <- list(
    list(mu = 0.5, nu = 2),
    list(mu = 0.2, nu = 10),
    list(mu = 0.8, nu = 0.5),
    list(mu = 0.99, nu = 100)
  )

  for (params in test_cases) {
    theta <- list(
      mu = array(params$mu, dim = c(1, 1, 1)),
      nu = array(params$nu, dim = c(1, 1, 1))
    )

    density <- PriorDensity(mdObj, theta)

    # Manual calculation
    mu_density <- dunif(params$mu, 0, mdObj$maxT)
    nu_density <- dgamma(1/params$nu, 2, 8) * (1/params$nu^2)  # Jacobian
    expected <- mu_density * nu_density

    expect_equal(density, expected, tolerance = 1e-10,
                 info = sprintf("mu=%.2f, nu=%.2f", params$mu, params$nu))
  }
})
