context("Weibull C++ and R implementations equivalence")

test_that("Weibull C++ and R implementations produce equivalent results", {
  # Generate test data
  set.seed(123)
  data <- rweibull(50, shape = 2, scale = 1)

  # Define prior parameters
  g0Priors <- c(10, 2, 4)  # phi, alpha0, beta0

  # Test with R backend
  set_use_cpp(FALSE)
  set.seed(456)
  dp_r <- DirichletProcessWeibull(data, g0Priors)
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)

  # Test with C++ backend (if available)
  skip_if(!can_use_cpp(DirichletProcessWeibull(data, g0Priors)),
          "C++ backend not available for Weibull")

  set_use_cpp(TRUE)
  set.seed(456)
  dp_cpp <- DirichletProcessWeibull(data, g0Priors)
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)

  # Compare results
  expect_equal(dp_r$numberClusters, dp_cpp$numberClusters, tolerance = 1)
  expect_equal(dp_r$alpha, dp_cpp$alpha, tolerance = 0.1)

  # Compare cluster parameters (allowing for minor numerical differences)
  if (dp_r$numberClusters == dp_cpp$numberClusters) {
    for (i in 1:2) {  # alpha and lambda parameters
      expect_equal(
        sort(as.numeric(dp_r$clusterParameters[[i]])),
        sort(as.numeric(dp_cpp$clusterParameters[[i]])),
        tolerance = 0.1
      )
    }
  }
})

test_that("Weibull C++ implementation provides significant speedup", {
  skip_if(!requireNamespace("microbenchmark", quietly = TRUE),
          "microbenchmark package not available")

  # Define prior parameters
  g0Priors <- c(10, 2, 4)

  # Check if C++ implementation is available
  skip_if(!can_use_cpp(DirichletProcessWeibull(rnorm(10), g0Priors)),
          "C++ implementation not available")

  # Generate larger dataset for benchmarking
  set.seed(789)
  data <- rweibull(200, shape = 2.5, scale = 1.5)

  # Benchmark both implementations
  bench_result <- microbenchmark::microbenchmark(
    R_implementation = {
      set_use_cpp(FALSE)
      dp <- DirichletProcessWeibull(data, g0Priors)
      Fit(dp, 50, progressBar = FALSE)
    },
    Cpp_implementation = {
      set_use_cpp(TRUE)
      dp <- DirichletProcessWeibull(data, g0Priors)
      Fit(dp, 50, progressBar = FALSE)
    },
    times = 5
  )

  # Calculate speedup
  mean_times <- aggregate(bench_result$time,
                          by = list(bench_result$expr),
                          FUN = mean)
  speedup <- mean_times$x[1] / mean_times$x[2]  # R time / C++ time

  # C++ should be significantly faster
  expect_gt(speedup, 2)  # At least 2x faster

  # Print speedup for information
  cat("\nWeibull C++ speedup:", round(speedup, 1), "x\n")
})

test_that("Weibull parameter estimation is accurate", {
  # Define prior parameters
  g0Priors <- c(10, 2, 4)

  # Check if C++ implementation is available
  skip_if(!can_use_cpp(DirichletProcessWeibull(rnorm(10), g0Priors)),
          "C++ implementation not available")

  # Generate data from known Weibull distribution
  set.seed(321)
  true_shape <- 2.0
  true_scale <- 1.5
  n_obs <- 100
  data <- rweibull(n_obs, shape = true_shape, scale = true_scale)

  # Fit with both implementations
  n_iter <- 200

  # R implementation
  set_use_cpp(FALSE)
  set.seed(654)
  dp_r <- DirichletProcessWeibull(data, g0Priors)
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  # C++ implementation
  set_use_cpp(TRUE)
  set.seed(654)
  dp_cpp <- DirichletProcessWeibull(data, g0Priors)
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  # Extract dominant cluster parameters (assuming most data in one cluster)
  get_dominant_params <- function(dp) {
    cluster_sizes <- table(dp$clusterLabels)
    dominant_cluster <- as.numeric(names(cluster_sizes)[which.max(cluster_sizes)])
    dominant_idx <- which(unique(dp$clusterLabels) == dominant_cluster)

    list(
      alpha = as.numeric(dp$clusterParameters[[1]][,,dominant_idx]),
      lambda = as.numeric(dp$clusterParameters[[2]][,,dominant_idx])
    )
  }

  params_r <- get_dominant_params(dp_r)
  params_cpp <- get_dominant_params(dp_cpp)

  # Convert lambda parameterization to scale
  # In the package: lambda = 1/scale^alpha
  scale_r <- (1/params_r$lambda)^(1/params_r$alpha)
  scale_cpp <- (1/params_cpp$lambda)^(1/params_cpp$alpha)

  # Check that estimates are reasonable (within 50% of true values)
  expect_lt(abs(params_r$alpha - true_shape) / true_shape, 0.5)
  expect_lt(abs(params_cpp$alpha - true_shape) / true_shape, 0.5)
  expect_lt(abs(scale_r - true_scale) / true_scale, 0.5)
  expect_lt(abs(scale_cpp - true_scale) / true_scale, 0.5)

  # R and C++ should produce similar estimates
  expect_equal(params_r$alpha, params_cpp$alpha, tolerance = 0.2)
  expect_equal(params_r$lambda, params_cpp$lambda, tolerance = 0.2)
})

test_that("Weibull handles edge cases correctly", {
  g0Priors <- c(10, 2, 4)

  # Test with very small dataset
  small_data <- rweibull(5, shape = 1, scale = 1)
  expect_error(
    dp_small <- DirichletProcessWeibull(small_data, g0Priors),
    NA  # Should not error
  )

  # Test with identical values
  identical_data <- rep(1.5, 20)
  expect_error(
    dp_identical <- DirichletProcessWeibull(identical_data, g0Priors),
    NA  # Should not error
  )

  # If C++ is available, test it too
  if (can_use_cpp(DirichletProcessWeibull(small_data, g0Priors))) {
    set_use_cpp(TRUE)

    expect_error(
      dp_small_cpp <- DirichletProcessWeibull(small_data, g0Priors),
      NA
    )

    expect_error(
      dp_identical_cpp <- DirichletProcessWeibull(identical_data, g0Priors),
      NA
    )
  }
})
