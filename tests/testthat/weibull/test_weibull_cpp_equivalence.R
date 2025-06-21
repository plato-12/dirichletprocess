context("Weibull Distribution C++ vs R Equivalence Tests")

test_that("Weibull C++ and R implementations produce equivalent results", {
  # Generate Weibull test data
  set.seed(42)
  true_alpha <- 2.0
  true_lambda <- 1.5
  n <- 100
  data <- rweibull(n, shape = true_alpha, scale = true_lambda)

  # Test R backend
  set_use_cpp(FALSE)
  set.seed(123)
  dp_r <- DirichletProcessWeibull(data)
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)

  # Test C++ backend (if available)
  skip_if(!can_use_cpp(DirichletProcessWeibull(data)),
          "C++ backend not available for Weibull")

  set_use_cpp(TRUE)
  set.seed(123)  # Same seed for reproducibility
  dp_cpp <- DirichletProcessWeibull(data)
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)

  # Compare results - allow for some stochastic variation
  cluster_diff <- abs(dp_r$numberClusters - dp_cpp$numberClusters)
  expect_true(cluster_diff <= 2,
              info = sprintf("Cluster difference: %d", cluster_diff))

  # Compare alpha values
  alpha_diff <- abs(dp_r$alpha - dp_cpp$alpha)
  expect_true(alpha_diff < 0.5,
              info = sprintf("Alpha difference: %.3f", alpha_diff))
})

test_that("Weibull C++ implementation provides significant speedup", {
  skip_if(!can_use_cpp(DirichletProcessWeibull(rnorm(10))),
          "C++ backend not available for Weibull")

  # Generate larger dataset for performance testing
  set.seed(42)
  data <- rweibull(500, shape = 2, scale = 1)

  # Time R implementation
  set_use_cpp(FALSE)
  r_time <- system.time({
    dp_r <- DirichletProcessWeibull(data)
    dp_r <- Fit(dp_r, 50, progressBar = FALSE)
  })["elapsed"]

  # Time C++ implementation
  set_use_cpp(TRUE)
  cpp_time <- system.time({
    dp_cpp <- DirichletProcessWeibull(data)
    dp_cpp <- Fit(dp_cpp, 50, progressBar = FALSE)
  })["elapsed"]

  speedup <- r_time / cpp_time

  expect_true(speedup > 10,
              info = sprintf("Speedup: %.1fx", speedup))
})

test_that("Weibull parameter estimation is accurate", {
  skip_if(!can_use_cpp(DirichletProcessWeibull(rnorm(10))),
          "C++ backend not available for Weibull")

  # Generate data from known Weibull distribution
  set.seed(42)
  true_alpha <- 3.0
  true_lambda <- 2.0
  data <- rweibull(200, shape = true_alpha, scale = true_lambda)

  # Fit using C++
  set_use_cpp(TRUE)
  dp <- DirichletProcessWeibull(data)
  dp <- Fit(dp, 200, progressBar = FALSE)

  # Extract cluster parameters
  # Assuming single cluster dominates
  main_cluster <- which.max(table(dp$clusterLabels))
  alpha_est <- dp$clusterParameters$alpha[main_cluster]
  lambda_est <- dp$clusterParameters$lambda[main_cluster]

  # Check parameter recovery (allowing for estimation error)
  expect_true(abs(alpha_est - true_alpha) < 0.5,
              info = sprintf("Alpha: true=%.2f, est=%.2f",
                             true_alpha, alpha_est))
  expect_true(abs(lambda_est - true_lambda) < 0.5,
              info = sprintf("Lambda: true=%.2f, est=%.2f",
                             true_lambda, lambda_est))
})
