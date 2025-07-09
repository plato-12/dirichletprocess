context("Exponential Distribution C++ Implementation Tests")

test_that("Exponential C++ likelihood matches R implementation", {
  x <- rexp(100, rate = 2)
  lambda <- 2

  # R implementation
  lik_r <- dexp(x, rate = lambda)

  # C++ implementation
  log_lik_cpp <- exponential_log_likelihood_cpp(x, lambda)
  lik_cpp <- exp(log_lik_cpp)

  expect_equal(lik_r, lik_cpp, tolerance = 1e-10)
})

test_that("Exponential posterior parameters are computed correctly", {
  # Generate test data
  set.seed(123)
  true_rate <- 3
  data <- matrix(rexp(50, rate = true_rate), ncol = 1)

  # Prior parameters
  alpha0 <- 2
  beta0 <- 1
  prior_params <- c(alpha0, beta0)

  # Expected posterior parameters
  n <- nrow(data)
  sum_x <- sum(data)
  expected_alpha <- alpha0 + n
  expected_beta <- beta0 + sum_x

  # C++ computation
  post_params <- exponential_posterior_parameters_cpp(prior_params, data)

  expect_equal(post_params$alpha, expected_alpha)
  expect_equal(post_params$beta, expected_beta)
})

test_that("Exponential C++ and R produce similar MCMC results", {
  # Generate data from mixture of exponentials
  set.seed(42)
  data <- c(rexp(30, rate = 2), rexp(30, rate = 5))

  # Test R backend
  set_use_cpp(FALSE)
  set.seed(123)
  dp_r <- DirichletProcessExponential(data)
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)

  # Test C++ backend
  skip_if(!can_use_cpp(DirichletProcessExponential(data)),
          "C++ backend not available for Exponential")

  set_use_cpp(TRUE)
  set.seed(123)
  dp_cpp <- DirichletProcessExponential(data)
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)

  # Compare number of clusters (should be close)
  cluster_diff <- abs(dp_r$numberClusters - dp_cpp$numberClusters)
  expect_true(cluster_diff <= 2,
              info = sprintf("R: %d clusters, C++: %d clusters",
                             dp_r$numberClusters, dp_cpp$numberClusters))

  # Both should find 2-3 clusters for this data
  expect_true(dp_r$numberClusters >= 1 && dp_r$numberClusters <= 4)
  expect_true(dp_cpp$numberClusters >= 1 && dp_cpp$numberClusters <= 4)
})

test_that("Exponential prior and posterior draws work correctly", {
  prior_params <- c(2, 1)  # alpha0 = 2, beta0 = 1

  # Test prior draw
  set.seed(123)
  prior_draw <- exponential_prior_draw_cpp(prior_params)
  # Extract the lambda value from the list
  lambda_value <- prior_draw$lambda[1,1,1]  # Get the single value from the 3D array
  expect_true(lambda_value > 0)
  expect_true(is.finite(lambda_value))

  # Test posterior draw
  data <- matrix(rexp(20, rate = 3), ncol = 1)
  post_draw <- exponential_posterior_draw_cpp(prior_params, data)
  # Extract the lambda value from the list
  post_lambda_value <- post_draw$lambda[1,1,1]  # Get the single value from the 3D array
  expect_true(post_lambda_value > 0)
  expect_true(is.finite(post_lambda_value))
})

test_that("Exponential handles edge cases properly", {
  # Test with very small rates
  x <- c(1000, 2000, 3000)  # Large values
  lambda <- 0.001  # Small rate

  log_lik <- exponential_log_likelihood_cpp(x, lambda)
  expect_true(all(is.finite(log_lik)))
  expect_true(all(log_lik < 0))  # Log likelihood should be negative

  # Test with zero/negative data (should return -Inf for negative)
  x_invalid <- c(-1, 0, 1)
  log_lik_invalid <- exponential_log_likelihood_cpp(x_invalid, 2)
  expect_equal(log_lik_invalid[1], -Inf)  # Negative value
  expect_true(is.finite(log_lik_invalid[2]))  # Zero is valid for exponential
  expect_true(is.finite(log_lik_invalid[3]))  # Positive value
})

test_that("Exponential C++ performance is faster than R", {
  skip_if(!requireNamespace("microbenchmark", quietly = TRUE))
  skip_if(!can_use_cpp(DirichletProcessExponential(rexp(100))),
          "C++ not available")

  data <- rexp(500, rate = 3)

  mb_result <- microbenchmark::microbenchmark(
    R_Backend = {
      set_use_cpp(FALSE)
      dp <- DirichletProcessExponential(data)
      Fit(dp, 50, progressBar = FALSE)
    },
    CPP_Backend = {
      set_use_cpp(TRUE)
      dp <- DirichletProcessExponential(data)
      Fit(dp, 50, progressBar = FALSE)
    },
    times = 3
  )

  # C++ should be significantly faster
  r_time <- median(mb_result$time[mb_result$expr == "R_Backend"])
  cpp_time <- median(mb_result$time[mb_result$expr == "CPP_Backend"])
  speedup <- r_time / cpp_time

  expect_true(speedup > 10,
              info = sprintf("C++ speedup: %.1fx", speedup))
})
