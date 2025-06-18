context("Beta C++ Implementation")

test_that("C++ Beta likelihood matches R implementation exactly", {
  skip_if_not(exists("beta_likelihood_cpp"), "C++ Beta not compiled")

  # Test various parameter combinations
  test_cases <- expand.grid(
    mu = c(0.1, 0.5, 0.9),
    tau = c(0.5, 5, 50),
    x = seq(0.05, 0.95, by = 0.1)
  )

  maxT <- 1

  for (i in 1:nrow(test_cases)) {
    mu <- test_cases$mu[i]
    tau <- test_cases$tau[i]
    x <- test_cases$x[i]

    # R calculation
    a <- (mu * tau) / maxT
    b <- (1 - mu/maxT) * tau
    expected <- (1/maxT) * dbeta(x/maxT, a, b)

    # C++ calculation
    actual <- beta_likelihood_cpp(x, mu, tau, maxT)

    expect_equal(actual, expected, tolerance = 1e-15,
                 info = sprintf("mu=%.2f, tau=%.2f, x=%.2f", mu, tau, x))
  }
})

test_that("C++ Beta handles numerical edge cases", {
  skip_if_not(exists("beta_likelihood_cpp"), "C++ Beta not compiled")

  # Test extreme tau values
  expect_true(is.finite(beta_likelihood_cpp(0.5, 0.5, 1e-10, 1)))
  expect_true(is.finite(beta_likelihood_cpp(0.5, 0.5, 1e10, 1)))

  # Test boundary x values
  expect_true(beta_likelihood_cpp(0, 0.5, 2, 1) >= 0)
  expect_true(beta_likelihood_cpp(1, 0.5, 2, 1) >= 0)

  # Test invalid parameters
  expect_equal(beta_likelihood_cpp(0.5, -0.1, 2, 1), 0)  # negative mu
  expect_equal(beta_likelihood_cpp(0.5, 1.1, 2, 1), 0)   # mu > maxT
  expect_equal(beta_likelihood_cpp(0.5, 0.5, -1, 1), 0) # negative tau
})

test_that("C++ Beta prior draw matches R distribution", {
  skip_if_not(exists("beta_prior_draw_cpp"), "C++ Beta not compiled")

  set.seed(123)
  n <- 10000
  priorParams <- c(2, 8)
  maxT <- 1

  # C++ draws
  cpp_draws <- beta_prior_draw_cpp(priorParams, maxT, n)
  mu_cpp <- as.vector(cpp_draws$mu)
  nu_cpp <- as.vector(cpp_draws$nu)

  # R draws
  set.seed(123)
  mu_r <- runif(n, 0, maxT)
  nu_r <- 1/rgamma(n, shape = priorParams[1], rate = priorParams[2])

  # Compare distributions using KS test
  ks_mu <- ks.test(mu_cpp, mu_r)
  ks_nu <- ks.test(nu_cpp, nu_r)

  expect_true(ks_mu$p.value > 0.01)
  expect_true(ks_nu$p.value > 0.01)
})
