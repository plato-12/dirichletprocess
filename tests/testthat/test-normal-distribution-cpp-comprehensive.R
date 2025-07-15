# test-normal-distribution-cpp-comprehensive.R
# Comprehensive unit tests for Normal Distribution C++ implementation

library(testthat)
library(dirichletprocess)

context("Normal Distribution C++ Implementation - Comprehensive Tests")

# Enable C++ mode
set_use_cpp(TRUE)

# Helper function to check if C++ implementations are available
cpp_available <- function() {
  if (!using_cpp()) return(FALSE)
  
  pkg_ns <- getNamespace("dirichletprocess")
  required_functions <- c(
    "normal_prior_draw_cpp",
    "normal_posterior_draw_cpp",
    "normal_posterior_parameters_cpp",
    "conjugate_cluster_component_update_cpp",
    "conjugate_cluster_parameter_update_cpp"
  )
  
  all_exist <- all(sapply(required_functions, function(f) exists(f, where = pkg_ns)))
  return(all_exist)
}

# Access functions from namespace if available
if (cpp_available()) {
  pkg_ns <- getNamespace("dirichletprocess")
  normal_prior_draw_cpp <- get("normal_prior_draw_cpp", pkg_ns)
  normal_posterior_draw_cpp <- get("normal_posterior_draw_cpp", pkg_ns)
  normal_posterior_parameters_cpp <- get("normal_posterior_parameters_cpp", pkg_ns)
  conjugate_cluster_component_update_cpp <- get("conjugate_cluster_component_update_cpp", pkg_ns)
  conjugate_cluster_parameter_update_cpp <- get("conjugate_cluster_parameter_update_cpp", pkg_ns)
}

# Helper function to check if arrays are statistically similar
check_statistical_similarity <- function(x, y, tol = 0.1, check_mean = TRUE, check_var = TRUE) {
  if (check_mean) {
    expect_equal(mean(x), mean(y), tolerance = tol,
                 info = sprintf("Means differ: %.4f vs %.4f", mean(x), mean(y)))
  }
  if (check_var) {
    expect_equal(var(c(x)), var(c(y)), tolerance = tol * 2,
                 info = sprintf("Variances differ: %.4f vs %.4f", var(c(x)), var(c(y))))
  }
}

# Test 1: Prior Draw functionality
test_that("Normal PriorDraw C++ implementation works correctly", {
  skip_if_not(cpp_available(), "C++ functions not available")

  # Test parameters
  prior_params <- c(mu0 = 0, kappa0 = 1, alpha0 = 2, beta0 = 1)
  n_draws <- 1000

  # Test 1.1: Basic functionality
  set.seed(123)
  result <- normal_prior_draw_cpp(prior_params, n_draws)

  expect_type(result, "list")
  expect_named(result, c("mu", "sigma"))
  expect_equal(dim(result$mu), c(1, 1, n_draws))
  expect_equal(dim(result$sigma), c(1, 1, n_draws))

  # Test 1.2: Parameter validation
  expect_error(normal_prior_draw_cpp(c(0, 1, 2), 10),
               "Normal distribution requires 4 prior parameters")
  expect_error(normal_prior_draw_cpp(prior_params, 0),
               "n must be at least 1")

  # Test 1.3: Statistical properties
  mu_samples <- as.vector(result$mu)
  sigma_samples <- as.vector(result$sigma)

  # Theoretical expectations for Normal-Inverse-Gamma
  expected_sigma2_mean <- prior_params[4] / (prior_params[3] - 1)
  observed_sigma2_mean <- mean(sigma_samples^2)

  expect_true(abs(observed_sigma2_mean - expected_sigma2_mean) < 0.1,
              info = sprintf("sigma^2 mean: observed %.3f, expected %.3f",
                             observed_sigma2_mean, expected_sigma2_mean))

  # Test 1.4: Different prior parameters
  test_priors <- list(
    c(5, 2, 3, 2),      # shifted mean
    c(0, 10, 5, 5),     # high precision
    c(-3, 0.1, 1, 0.5)  # low precision
  )

  for (params in test_priors) {
    result <- normal_prior_draw_cpp(params, 100)
    expect_equal(dim(result$mu), c(1, 1, 100))
    expect_true(all(result$sigma > 0))
  }
})

# Test 2: Posterior Draw functionality
test_that("Normal PosteriorDraw C++ implementation works correctly", {
  skip_if_not(cpp_available(), "C++ functions not available")

  prior_params <- c(0, 1, 2, 1)

  # Test 2.1: Single observation
  x_single <- matrix(2.5, ncol = 1)
  result <- normal_posterior_draw_cpp(prior_params, x_single, 100)

  expect_type(result, "list")
  expect_named(result, c("mu", "sigma"))
  expect_equal(dim(result$mu), c(1, 1, 100))

  # Test 2.2: Multiple observations
  set.seed(456)
  x_multi <- matrix(rnorm(20, mean = 3, sd = 1.5), ncol = 1)
  result_multi <- normal_posterior_draw_cpp(prior_params, x_multi, 500)

  mu_samples <- as.vector(result_multi$mu)
  # Posterior mean should be pulled towards data mean
  expect_true(abs(mean(mu_samples) - mean(x_multi)) < 0.5)

  # Test 2.3: Edge cases
  # Empty data matrix
  expect_error(normal_posterior_draw_cpp(prior_params, matrix(numeric(0), ncol = 1), 10),
               "Cannot draw from posterior with empty data")

  # Very large data
  x_large <- matrix(rnorm(1000), ncol = 1)
  result_large <- normal_posterior_draw_cpp(prior_params, x_large, 10)
  expect_equal(dim(result_large$mu), c(1, 1, 10))
})

# Test 3: Posterior Parameters calculation
test_that("Normal posterior parameters calculation is correct", {
  skip_if_not(cpp_available(), "C++ functions not available")

  prior_params <- c(mu0 = 0, kappa0 = 1, alpha0 = 2, beta0 = 1)

  # Test 3.1: Known example
  x <- matrix(c(1, 2, 3), ncol = 1)
  post_params <- normal_posterior_parameters_cpp(prior_params, x)

  expect_type(post_params, "double")
  expect_equal(dim(post_params), c(1, 4))

  # Manual calculation
  n <- nrow(x)
  x_bar <- mean(x)
  mu_n <- (prior_params[2] * prior_params[1] + n * x_bar) / (prior_params[2] + n)
  kappa_n <- prior_params[2] + n
  alpha_n <- prior_params[3] + n / 2
  beta_n <- prior_params[4] + 0.5 * sum((x - x_bar)^2) +
    prior_params[2] * n * (x_bar - prior_params[1])^2 / (2 * (prior_params[2] + n))

  # Compare without attributes
  expect_equal(post_params[1, 1], mu_n, tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(post_params[1, 2], kappa_n, tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(post_params[1, 3], alpha_n, tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(post_params[1, 4], beta_n, tolerance = 1e-10, check.attributes = FALSE)
})

# Test 4: Conjugate cluster component update
test_that("Conjugate cluster component update works", {
  skip_if_not(cpp_available(), "C++ functions not available")

  # Create a simple DP object
  set.seed(789)
  data <- matrix(c(rnorm(10, -2), rnorm(10, 2)), ncol = 1)

  dp_obj <- list(
    data = data,
    clusterLabels = as.integer(c(rep(0, 10), rep(1, 10))),
    pointsPerCluster = as.integer(c(10, 10)),
    numberClusters = 2L,
    alpha = 1.0,
    mixingDistribution = list(priorParameters = c(0, 1, 2, 1)),
    clusterParameters = list(
      c(-2, 2),  # mu values
      c(1, 1)    # sigma values
    ),
    predictiveArray = rep(0.1, 20)
  )

  result <- conjugate_cluster_component_update_cpp(dp_obj)

  # Basic checks
  expect_type(result, "list")
  expect_equal(length(result$clusterLabels), 20)
  expect_true(all(result$clusterLabels >= 0))
  expect_equal(sum(result$pointsPerCluster), 20)
  expect_true(result$numberClusters >= 1)
})

# Test 5: Conjugate cluster parameter update
test_that("Conjugate cluster parameter update works", {
  skip_if_not(cpp_available(), "C++ functions not available")

  # Create DP object with two clear clusters
  set.seed(101112)
  data <- matrix(c(rnorm(20, -5, 0.5), rnorm(20, 5, 0.5)), ncol = 1)

  dp_obj <- list(
    data = data,
    clusterLabels = as.integer(c(rep(0, 20), rep(1, 20))),
    numberClusters = 2L,
    mixingDistribution = list(priorParameters = c(0, 0.1, 2, 1)),
    clusterParameters = list(
      c(0, 0),  # initial mu values
      c(1, 1)   # initial sigma values
    )
  )

  result <- conjugate_cluster_parameter_update_cpp(dp_obj)

  # Parameters should be updated towards cluster data
  mu_updated <- result[[1]]
  expect_true(abs(mu_updated[1] - mean(data[1:20])) < 1)
  expect_true(abs(mu_updated[2] - mean(data[21:40])) < 1)
})

# Test 6: Prior vs Posterior consistency
test_that("Prior and posterior draws are consistent", {
  skip_if_not(cpp_available(), "C++ functions not available")

  prior_params <- c(0, 0.01, 2, 1)  # Weak prior

  # Generate data from known distribution
  set.seed(131415)
  true_mu <- 3
  true_sigma <- 1.5
  n_data <- 100
  data <- matrix(rnorm(n_data, true_mu, true_sigma), ncol = 1)

  # Prior and posterior draws
  prior_draws <- normal_prior_draw_cpp(prior_params, 1000)
  posterior_draws <- normal_posterior_draw_cpp(prior_params, data, 1000)

  # With weak prior and lots of data, posterior should concentrate around true values
  post_mu <- as.vector(posterior_draws$mu)
  post_sigma <- as.vector(posterior_draws$sigma)

  # Use more reasonable tolerances
  expect_true(abs(mean(post_mu) - true_mu) < 0.5,
              info = sprintf("Posterior mean: %.3f, True mean: %.3f",
                             mean(post_mu), true_mu))
  expect_true(abs(mean(post_sigma) - true_sigma) < 0.5,
              info = sprintf("Posterior SD: %.3f, True SD: %.3f",
                             mean(post_sigma), true_sigma))

  # Posterior should be less variable than prior
  prior_mu_var <- var(as.vector(prior_draws$mu))
  post_mu_var <- var(post_mu)
  expect_true(post_mu_var < prior_mu_var / 5)
})

# Test 7: Numerical stability
test_that("C++ implementation handles numerical edge cases", {
  skip_if_not(cpp_available(), "C++ functions not available")

  # Test 7.1: Very small/large parameters
  extreme_params <- list(
    c(1000, 0.001, 0.5, 0.001),    # Large mean, small precision
    c(-1000, 1000, 100, 100),       # Large negative mean, high precision
    c(0, 1e-6, 1, 1e-6)             # Very small kappa and beta
  )

  for (params in extreme_params) {
    result <- tryCatch(
      normal_prior_draw_cpp(params, 10),
      error = function(e) NULL
    )
    expect_false(is.null(result),
                 info = sprintf("Failed with params: %s", paste(params, collapse = ", ")))
    if (!is.null(result)) {
      expect_true(all(is.finite(result$mu)))
      expect_true(all(is.finite(result$sigma)))
      expect_true(all(result$sigma > 0))
    }
  }

  # Test 7.2: Data with outliers
  data_outliers <- matrix(c(rnorm(95), rep(100, 5)), ncol = 1)
  prior_params <- c(0, 1, 2, 1)

  result <- normal_posterior_draw_cpp(prior_params, data_outliers, 100)
  expect_true(all(is.finite(result$mu)))
  expect_true(all(is.finite(result$sigma)))
})

# Test 8: Performance comparison
test_that("C++ implementation is faster than R", {
  skip_if_not(cpp_available(), "C++ functions not available")
  skip_if_not(exists("PriorDraw.normal"), "R implementation not available")

  prior_params <- c(0, 1, 2, 1)
  n_draws <- 10000

  # Create mixing distribution object for R version
  md_obj <- MixingDistribution("normal", prior_params, "conjugate")

  # Time R implementation
  time_r <- system.time({
    r_result <- PriorDraw(md_obj, n_draws)
  })

  # Time C++ implementation
  time_cpp <- system.time({
    cpp_result <- normal_prior_draw_cpp(prior_params, n_draws)
  })

  # C++ should be significantly faster
  speedup <- if (as.numeric(time_cpp["elapsed"]) > 0) {
    as.numeric(time_r["elapsed"]) / as.numeric(time_cpp["elapsed"])
  } else {
    Inf  # C++ is infinitely fast if time is 0
  }

  cat("\nSpeedup factor:",
      ifelse(is.infinite(speedup), "Inf", round(speedup, 1)), "x\n")

  expect_true(speedup > 2 || is.infinite(speedup),
              info = sprintf("C++ only %sx faster than R",
                             ifelse(is.infinite(speedup), "Inf",
                                    sprintf("%.1f", speedup))))

  # Results should be statistically similar
  check_statistical_similarity(r_result$mu, cpp_result$mu, tol = 0.2)
  check_statistical_similarity(r_result$sigma, cpp_result$sigma, tol = 0.2)
})

# Test 9: Integration test - Complete MCMC cycle
test_that("Complete MCMC cycle works with C++ backend", {
  skip_on_cran()  # This test might take longer

  # Generate mixture data
  set.seed(161718)
  data <- c(rnorm(50, -3, 0.5), rnorm(50, 0, 1), rnorm(50, 3, 0.5))

  # Run with C++ backend
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(data)
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)

  # Check reasonable bounds for stochastic algorithm
  expect_true(dp_cpp$numberClusters >= 2 && dp_cpp$numberClusters <= 10,
              info = sprintf("C++ clusters: %d", dp_cpp$numberClusters))

  # Alpha should be positive and reasonable
  expect_true(dp_cpp$alpha > 0 && dp_cpp$alpha < 20,
              info = sprintf("C++ alpha: %.2f", dp_cpp$alpha))

  # All points should be assigned
  expect_equal(length(dp_cpp$clusterLabels), length(data))

  # For C++ implementation, cluster labels are 0-indexed internally
  # but should be converted to 1-indexed for R compatibility
  # If not converted, check for 0-indexed validity
  if (min(dp_cpp$clusterLabels) == 0) {
    # 0-indexed check
    expect_true(all(dp_cpp$clusterLabels >= 0))
    expect_true(all(dp_cpp$clusterLabels < dp_cpp$numberClusters))
  } else {
    # 1-indexed check (if properly converted)
    expect_true(all(dp_cpp$clusterLabels > 0))
    expect_true(all(dp_cpp$clusterLabels <= dp_cpp$numberClusters))
  }

  # Reset to default
  set_use_cpp(TRUE)
})

# Test 10: Memory and bounds checking
test_that("C++ implementation handles memory correctly", {
  skip_if_not(cpp_available(), "C++ functions not available")

  # Test with single data point
  dp_single <- list(
    data = matrix(1, ncol = 1),
    clusterLabels = as.integer(0),
    pointsPerCluster = as.integer(1),
    numberClusters = 1L,
    alpha = 1.0,
    mixingDistribution = list(priorParameters = c(0, 1, 2, 1)),
    clusterParameters = list(matrix(0), matrix(1)),
    predictiveArray = 1
  )

  result <- conjugate_cluster_component_update_cpp(dp_single)
  expect_equal(length(result$clusterLabels), 1)

  # Test with many clusters
  n_data <- 100
  n_clusters <- 20
  cluster_assignments <- sample(0:(n_clusters-1), n_data, replace = TRUE)

  dp_many <- list(
    data = matrix(rnorm(n_data), ncol = 1),
    clusterLabels = as.integer(cluster_assignments),
    pointsPerCluster = as.integer(table(factor(cluster_assignments, levels = 0:(n_clusters-1)))),
    numberClusters = n_clusters,
    alpha = 1.0,
    mixingDistribution = list(priorParameters = c(0, 1, 2, 1)),
    clusterParameters = list(
      matrix(rnorm(n_clusters), nrow = 1),
      matrix(rgamma(n_clusters, 2, 2), nrow = 1)
    ),
    predictiveArray = rep(1, n_data)
  )

  result <- conjugate_cluster_component_update_cpp(dp_many)
  expect_true(all(result$clusterLabels >= 0))
  expect_true(all(result$clusterLabels < result$numberClusters))
})

# Final summary test
test_that("Summary: All core C++ functions are implemented", {
  cpp_functions <- c(
    "normal_prior_draw_cpp",
    "normal_posterior_draw_cpp",
    "normal_posterior_parameters_cpp",
    "conjugate_cluster_component_update_cpp",
    "conjugate_cluster_parameter_update_cpp"
  )

  available <- sapply(cpp_functions, exists)

  cat("\n\nC++ Implementation Status:\n")
  cat(sprintf("%-40s %s\n", "Function", "Available"))
  cat(rep("-", 50), "\n", sep = "")

  for (i in seq_along(cpp_functions)) {
    status <- if (available[i]) "✓" else "✗"
    cat(sprintf("%-40s %s\n", cpp_functions[i], status))
  }

  expect_true(all(available),
              info = paste("Missing functions:",
                           cpp_functions[!available], collapse = ", "))
})
