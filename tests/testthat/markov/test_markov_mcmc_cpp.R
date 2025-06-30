# tests/testthat/test_markov_mcmc_cpp.R
context("Markov MCMC C++ Implementation")

test_that("Markov MCMC C++ produces valid results", {
  skip_if_not(exists("run_markov_mcmc_cpp"))

  # Generate Markov chain data
  set.seed(123)
  true_states <- c(rep(1, 50), rep(2, 50), rep(3, 50))
  data <- numeric(150)
  data[true_states == 1] <- rnorm(50, 0, 1)
  data[true_states == 2] <- rnorm(50, 5, 1)
  data[true_states == 3] <- rnorm(50, 10, 1)

  # Create DP object
  md <- GaussianMixtureCreate()
  dp <- DirichletHMMCreate(data, md, alpha = 1, beta = 1)

  # Run C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- run_markov_mcmc_cpp(dp, its = 100, update_prior = TRUE)

  # Check results
  expect_is(dp_cpp$states, "integer")
  expect_length(dp_cpp$states, 150)
  expect_true(all(dp_cpp$states >= 1))
  expect_true(length(unique(dp_cpp$states)) >= 2)  # At least 2 states
  expect_true(length(unique(dp_cpp$states)) <= 10) # Not too many states

  # Check chains
  expect_length(dp_cpp$alphaChain, 90)  # 100 iterations - 10% burn-in
  expect_length(dp_cpp$betaChain, 90)
  expect_true(all(dp_cpp$alphaChain > 0))
  expect_true(all(dp_cpp$betaChain > 0))
})

test_that("Markov MCMC handles edge cases", {
  skip_if_not(exists("run_markov_mcmc_cpp"))

  # Single observation
  data <- matrix(rnorm(1), ncol = 1)
  md <- GaussianMixtureCreate()
  dp <- DirichletHMMCreate(data, md, alpha = 1, beta = 1)

  expect_error(run_markov_mcmc_cpp(dp, its = 10), NA)

  # Two observations
  data <- matrix(rnorm(2), ncol = 1)
  dp <- DirichletHMMCreate(data, md, alpha = 1, beta = 1)
  result <- run_markov_mcmc_cpp(dp, its = 10)

  expect_length(result$states, 2)
})
