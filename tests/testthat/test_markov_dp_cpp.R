context("Markov DP C++ Implementation")

# Ensure the required functions are loaded
library(dirichletprocess)

# Helper function to create test data
create_hmm_test_data <- function(n = 150, states = c(1, 3, 5)) {
  # Generate data with clear state transitions
  data1 <- rnorm(n/3, states[1], sqrt(3))
  data2 <- rnorm(n/3, states[2], sqrt(3))
  data3 <- rnorm(n/3, states[3], sqrt(3))
  c(data1, data2, data3)
}

# Helper to create a basic Markov DP object
create_test_markov_dp <- function(data = NULL) {
  if (is.null(data)) {
    data <- create_hmm_test_data()
  }

  # Use the correct function name from the package
  mdObj <- GaussianMixtureCreate()
  dp <- DirichletHMMCreate(data, mdObj, alpha = 2, beta = 3)
  return(dp)
}

test_that("Markov DP can be created from R object", {
  dp <- create_test_markov_dp()

  # Test basic structure
  expect_is(dp, c("list", "markov", "dirichletprocess", "normal", "conjugate"))
  expect_equal(length(dp$data), 150)
  expect_equal(dp$alpha, 2)
  expect_equal(dp$beta, 3)
  expect_equal(length(dp$states), length(dp$data))
  expect_equal(length(dp$params), length(dp$data))
})

test_that("C++ Markov DP creation and conversion", {
  skip_if_not(exists("markov_dp_create_cpp"))

  dp_r <- create_test_markov_dp()

  # Create C++ object and convert back
  dp_cpp <- markov_dp_create_cpp(dp_r)

  expect_is(dp_cpp, "list")
  expect_true(!is.null(attr(dp_cpp, "cpp_ptr")))
})

test_that("UpdateStates C++ implementation works correctly", {
  skip_if_not(exists("markov_dp_update_states_cpp"))

  dp <- create_test_markov_dp()
  initial_states <- dp$states

  # Convert to 0-indexed for C++
  dp$states <- dp$states - 1

  # Update states using C++
  result <- markov_dp_update_states_cpp(dp)

  # Convert back to 1-indexed
  result$states <- result$states + 1

  # Check structure
  expect_is(result$states, "integer")
  expect_equal(length(result$states), length(dp$data))
  expect_true(all(result$states >= 1))

  # States should have changed (with high probability)
  expect_false(all(result$states == initial_states))

  # Check that params are updated correctly
  expect_is(result$params, "list")
  expect_equal(length(result$params), length(result$states))
})

test_that("UpdateAlphaBeta C++ implementation works correctly", {
  skip_if_not(exists("markov_dp_update_alpha_beta_cpp"))

  dp <- create_test_markov_dp()
  initial_alpha <- dp$alpha
  initial_beta <- dp$beta

  # Convert states to 0-indexed
  dp$states <- dp$states - 1

  # Update alpha and beta using C++
  result <- markov_dp_update_alpha_beta_cpp(dp)

  # Check that alpha and beta were updated
  expect_is(result$alpha, "numeric")
  expect_is(result$beta, "numeric")
  expect_true(result$alpha > 0)
  expect_true(result$beta > 0)

  # They should be different from initial values (with high probability)
  expect_false(result$alpha == initial_alpha && result$beta == initial_beta)
})

test_that("ParamUpdate C++ implementation works correctly", {
  skip_if_not(exists("markov_dp_param_update_cpp"))

  dp <- create_test_markov_dp()

  # Convert states to 0-indexed
  dp$states <- dp$states - 1

  # Update parameters using C++
  result <- markov_dp_param_update_cpp(dp)

  # Check structure
  expect_is(result$uniqueParams, "list")
  expect_is(result$params, "list")
  expect_equal(length(result$params), length(result$states))

  # Check that unique parameters match the number of unique states
  unique_states <- length(unique(result$states))
  if (length(result$uniqueParams) > 0) {
    first_param <- result$uniqueParams[[1]]
    if (length(dim(first_param)) == 3) {
      expect_equal(dim(first_param)[3], unique_states)
    }
  }
})

test_that("Full Markov DP fit works with C++", {
  skip_if_not(exists("markov_dp_fit_cpp"))

  dp <- create_test_markov_dp()

  # Convert states to 0-indexed
  dp$states <- dp$states - 1

  # Fit using C++
  iterations <- 10
  result <- markov_dp_fit_cpp(dp, iterations, updatePrior = FALSE, progressBar = FALSE)

  # Convert states back to 1-indexed
  result$states <- result$states + 1
  if (!is.null(result$statesChain)) {
    result$statesChain <- lapply(result$statesChain, function(x) x + 1)
  }

  # Check chains
  expect_equal(length(result$alphaChain), iterations)
  expect_equal(length(result$betaChain), iterations)
  expect_equal(length(result$statesChain), iterations)
  expect_equal(length(result$paramChain), iterations)

  # Check that values change over iterations
  expect_false(all(result$alphaChain == result$alphaChain[1]))
  expect_false(all(result$betaChain == result$betaChain[1]))

  # All chain values should be positive
  expect_true(all(result$alphaChain > 0))
  expect_true(all(result$betaChain > 0))
})

test_that("C++ and R implementations produce similar results", {
  skip_if_not(exists("markov_dp_fit_cpp"))
  skip_if_not(exists("fit_hmm"))

  set.seed(123)
  data <- create_hmm_test_data(n = 60)  # Smaller dataset for faster comparison

  # Create two identical DP objects
  dp_r <- create_test_markov_dp(data)
  dp_cpp <- create_test_markov_dp(data)

  # Set same initial states
  set.seed(456)
  init_states <- sample(1:3, length(data), replace = TRUE)
  dp_r$states <- init_states
  dp_cpp$states <- init_states

  # Fit with R implementation
  set.seed(789)
  dp_r_fit <- fit_hmm(dp_r, its = 5, progressBar = FALSE)

  # Fit with C++ implementation
  set.seed(789)
  dp_cpp$states <- dp_cpp$states - 1  # Convert to 0-indexed
  dp_cpp_fit <- markov_dp_fit_cpp(dp_cpp, iterations = 5, updatePrior = FALSE, progressBar = FALSE)
  dp_cpp_fit$states <- dp_cpp_fit$states + 1  # Convert back

  # Results should be similar (not identical due to implementation differences)
  # Check that final alpha/beta are in similar ranges
  expect_true(abs(dp_r_fit$alpha - dp_cpp_fit$alpha) / dp_r_fit$alpha < 0.5)
  expect_true(abs(dp_r_fit$beta - dp_cpp_fit$beta) / dp_r_fit$beta < 0.5)

  # Check that the number of unique states is similar
  n_states_r <- length(unique(dp_r_fit$states))
  n_states_cpp <- length(unique(dp_cpp_fit$states))
  expect_true(abs(n_states_r - n_states_cpp) <= 2)
})

test_that("C++ implementation handles edge cases", {
  skip_if_not(exists("markov_dp_update_states_cpp"))

  # Test with very small dataset
  small_data <- rnorm(5)
  dp_small <- DirichletHMMCreate(small_data, GaussianMixtureCreate(), 1, 1)
  dp_small$states <- dp_small$states - 1

  result_small <- markov_dp_update_states_cpp(dp_small)
  expect_equal(length(result_small$states), 5)
  expect_true(all(result_small$states >= 0))

  # Test with single state
  dp_single <- create_test_markov_dp()
  dp_single$states <- rep(0, length(dp_single$data))  # All same state (0-indexed)

  result_single <- markov_dp_update_states_cpp(dp_single)
  expect_equal(length(result_single$states), length(dp_single$data))

  # Test with maximum number of states
  dp_max <- create_test_markov_dp()
  dp_max$states <- seq(0, length(dp_max$data) - 1)  # Each point is its own state

  result_max <- markov_dp_update_states_cpp(dp_max)
  expect_equal(length(result_max$states), length(dp_max$data))
  expect_true(length(unique(result_max$states)) <= length(dp_max$data))
})

test_that("C++ wrapper functions work correctly", {
  skip_if_not(exists("using_cpp_markov_samplers"))

  # Test enabling/disabling C++ samplers
  old_setting <- using_cpp_markov_samplers()

  enable_cpp_markov_samplers(TRUE)
  expect_true(using_cpp_markov_samplers())

  enable_cpp_markov_samplers(FALSE)
  expect_false(using_cpp_markov_samplers())

  # Restore original setting
  enable_cpp_markov_samplers(old_setting)
})

test_that("Fit.markov dispatches correctly to C++", {
  skip_if_not(exists("Fit.markov.cpp"))

  dp <- create_test_markov_dp()

  # Enable C++ samplers
  old_setting <- using_cpp_markov_samplers()
  enable_cpp_markov_samplers(TRUE)

  # This should use C++ implementation
  result <- Fit(dp, its = 5)

  expect_is(result, "markov")
  expect_equal(length(result$alphaChain), 5)
  expect_equal(length(result$betaChain), 5)

  # Restore setting
  enable_cpp_markov_samplers(old_setting)
})

test_that("Parameter structure is maintained correctly", {
  skip_if_not(exists("markov_dp_param_update_cpp"))

  dp <- create_test_markov_dp()

  # Ensure states are diverse
  dp$states <- rep(1:3, length.out = length(dp$data))
  dp$states <- dp$states - 1  # Convert to 0-indexed

  result <- markov_dp_param_update_cpp(dp)

  # Check parameter structure for normal distribution
  expect_true(length(result$uniqueParams) >= 2)  # mu and sigma

  # Each unique state should have parameters
  unique_states <- unique(result$states)
  n_unique <- length(unique_states)

  for (param in result$uniqueParams) {
    if (length(dim(param)) == 3) {
      expect_equal(dim(param)[3], n_unique)
    }
  }

  # Check that params list references correct unique parameters
  for (i in seq_along(result$params)) {
    state_params <- result$params[[i]]
    expect_is(state_params, "list")
    expect_true(length(state_params) >= 2)  # mu and sigma for normal
  }
})

test_that("State relabeling works correctly", {
  skip_if_not(exists("markov_dp_update_states_cpp"))

  dp <- create_test_markov_dp()

  # Create non-contiguous states
  dp$states <- c(rep(0, 50), rep(2, 50), rep(5, 50))  # 0, 2, 5 (0-indexed)

  result <- markov_dp_update_states_cpp(dp)

  # States should be relabeled to be contiguous
  unique_states <- sort(unique(result$states))
  expect_equal(unique_states, seq(0, length(unique_states) - 1))
})
