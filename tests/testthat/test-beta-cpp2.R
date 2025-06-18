context("Beta Distribution C++ Implementation Tests")

test_that("Beta C++ likelihood matches R implementation", {
  # Test data
  x <- c(0.2, 0.5, 0.8)
  mu <- 0.5
  tau <- 5
  maxT <- 1

  # R implementation
  a <- (mu * tau) / maxT
  b <- (1 - mu/maxT) * tau
  expected <- (1/maxT) * dbeta(x/maxT, a, b)

  # C++ implementation (if available)
  skip_if_not(exists("_dirichletprocess_beta_likelihood_cpp"),
              "C++ Beta implementation not compiled")

  actual <- beta_likelihood_cpp(x, mu, tau, maxT)

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("Beta C++ and R produce equivalent results", {
  skip_if_not(can_use_cpp(DirichletProcessBeta(rbeta(10, 2, 2))),
              "C++ backend not available for Beta")

  # Generate test data
  set.seed(123)
  data <- c(rbeta(30, 2, 8), rbeta(30, 8, 2))

  # R backend
  set_use_cpp(FALSE)
  set.seed(456)
  dp_r <- DirichletProcessBeta(data, g0Priors = c(2, 8))
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)

  # C++ backend
  set_use_cpp(TRUE)
  set.seed(456)
  dp_cpp <- DirichletProcessBeta(data, g0Priors = c(2, 8))
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)

  # Compare results - allow for some MCMC variation
  expect_equal(dp_r$numberClusters, dp_cpp$numberClusters, tolerance = 2)
  expect_equal(mean(dp_r$clusterLabels), mean(dp_cpp$clusterLabels), tolerance = 0.5)
})

test_that("Beta handles edge cases correctly", {
  skip_if_not(can_use_cpp(DirichletProcessBeta(rbeta(10, 2, 2))),
              "C++ backend not available for Beta")

  # Extreme data near boundaries
  extreme_data <- c(rep(0.01, 10), rep(0.99, 10))

  set_use_cpp(TRUE)
  expect_error({
    dp <- DirichletProcessBeta(extreme_data)
    dp <- Fit(dp, 50, progressBar = FALSE)
  }, NA)  # Should not error

  expect_true(dp$numberClusters >= 1)
})

test_that("Beta C++ performance is faster than R", {
  skip_if_not(can_use_cpp(DirichletProcessBeta(rbeta(10, 2, 2))),
              "C++ backend not available for Beta")
  skip_if_not(requireNamespace("microbenchmark", quietly = TRUE))

  data <- rbeta(100, 3, 3)

  mb_result <- microbenchmark::microbenchmark(
    R_Backend = {
      set_use_cpp(FALSE)
      dp <- DirichletProcessBeta(data)
      Fit(dp, 50, progressBar = FALSE)
    },
    CPP_Backend = {
      set_use_cpp(TRUE)
      dp <- DirichletProcessBeta(data)
      Fit(dp, 50, progressBar = FALSE)
    },
    times = 5
  )

  r_time <- median(mb_result$time[mb_result$expr == "R_Backend"])
  cpp_time <- median(mb_result$time[mb_result$expr == "CPP_Backend"])

  # Expect at least 10x speedup
  expect_true(r_time / cpp_time > 10)
})
