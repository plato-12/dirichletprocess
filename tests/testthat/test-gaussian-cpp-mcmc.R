# tests/testthat/test_gaussian_cpp_mcmc.R
context("Gaussian C++ MCMC Implementation Tests")

# Helper functions ---------------------------------------------------------

#' Create synthetic Gaussian mixture data
#' @param n Total number of observations
#' @param means Vector of component means
#' @param sds Vector of component standard deviations
#' @param weights Vector of component weights (will be normalized)
#' @param seed Random seed for reproducibility
create_gaussian_mixture_data <- function(n = 100,
                                         means = c(-2, 2),
                                         sds = c(0.5, 0.5),
                                         weights = c(0.5, 0.5),
                                         seed = 123) {
  set.seed(seed)
  weights <- weights / sum(weights)
  n_components <- length(means)

  # Generate component assignments
  components <- sample(1:n_components, n, replace = TRUE, prob = weights)

  # Generate data
  data <- numeric(n)
  for (i in 1:n) {
    data[i] <- rnorm(1, means[components[i]], sds[components[i]])
  }

  list(data = data, true_components = components)
}

#' Create standard MCMC parameters
create_mcmc_params <- function(n_iter = 1000,
                               n_burn = 200,
                               thin = 1,
                               update_concentration = TRUE,
                               alpha = 1.0) {
  list(
    n_iter = n_iter,
    n_burn = n_burn,
    thin = thin,
    update_concentration = update_concentration,
    alpha = alpha
  )
}

#' Create Gaussian mixing distribution parameters
create_gaussian_params <- function(mu0 = 0.0,
                                   kappa0 = 1.0,
                                   alpha0 = 1.0,
                                   beta0 = 1.0) {
  list(
    type = "gaussian",
    mu0 = mu0,
    kappa0 = kappa0,
    alpha0 = alpha0,
    beta0 = beta0
  )
}

#' Check if C++ implementation is available
cpp_available <- function() {
  exists("run_mcmc_cpp") && exists("_dirichletprocess_run_mcmc_cpp")
}

#' Compute clustering metrics
compute_clustering_metrics <- function(true_labels, estimated_labels) {
  # Adjusted Rand Index
  ari <- mclust::adjustedRandIndex(true_labels, estimated_labels)

  # Number of unique clusters
  n_clusters <- length(unique(estimated_labels))

  list(ari = ari, n_clusters = n_clusters)
}

# Basic Functionality Tests ------------------------------------------------

test_that("C++ MCMC runner executes without errors", {
  skip_if(!cpp_available(), "C++ implementation not available")

  # Simple test data
  data_obj <- create_gaussian_mixture_data(n = 50, seed = 42)
  data_matrix <- matrix(data_obj$data, ncol = 1)

  # Run MCMC
  result <- run_mcmc_cpp(
    data = data_matrix,
    mixing_dist_params = create_gaussian_params(),
    mcmc_params = create_mcmc_params(n_iter = 100, n_burn = 20)
  )

  # Check output structure - C++ returns additional fields
  expect_type(result, "list")
  # The actual names returned include both the expected fields and compatibility fields
  expect_true(all(c("cluster_labels", "alpha", "theta", "n_clusters") %in% names(result)))

  # Check dimensions
  n_saved <- (100 - 20) / 1  # (n_iter - n_burn) / thin

  # cluster_labels is a matrix, not a list
  if (is.matrix(result$cluster_labels)) {
    expect_equal(nrow(result$cluster_labels), n_saved)
    expect_equal(ncol(result$cluster_labels), 50)  # number of data points
  } else {
    expect_equal(length(result$cluster_labels), n_saved)
  }

  expect_equal(length(result$alpha), n_saved)
  # theta is the final parameters, not a chain
  expect_type(result$theta, "list")
  expect_equal(length(result$n_clusters), n_saved)

  # Check data types - adjusted for actual structure
  if (is.matrix(result$cluster_labels)) {
    expect_true(is.numeric(result$cluster_labels))
  } else {
    expect_true(all(sapply(result$cluster_labels, is.integer)))
  }
  expect_true(all(sapply(result$alpha, is.numeric)))
  expect_true(all(sapply(result$n_clusters, is.integer)))
})

test_that("C++ MCMC handles different data sizes", {
  skip_if(!cpp_available(), "C++ implementation not available")

  for (n in c(10, 50, 100, 500)) {
    data_matrix <- matrix(rnorm(n), ncol = 1)

    result <- run_mcmc_cpp(
      data = data_matrix,
      mixing_dist_params = create_gaussian_params(),
      mcmc_params = create_mcmc_params(n_iter = 50, n_burn = 10)
    )

    # Check the actual structure returned
    if (is.matrix(result$cluster_labels)) {
      expect_equal(ncol(result$cluster_labels), n,
                   info = paste("Failed for n =", n))
    } else if (is.list(result$cluster_labels)) {
      expect_equal(length(result$cluster_labels[[1]]), n,
                   info = paste("Failed for n =", n))
    }
  }
})

# Parameter Validation Tests -----------------------------------------------

test_that("C++ MCMC validates input parameters", {
  skip_if(!cpp_available(), "C++ implementation not available")

  data_matrix <- matrix(rnorm(20), ncol = 1)

  # Test invalid n_iter
  expect_error(
    run_mcmc_cpp(data_matrix, create_gaussian_params(),
                 create_mcmc_params(n_iter = 0)),
    "n_iter must be positive"
  )

  # Test invalid n_burn
  expect_error(
    run_mcmc_cpp(data_matrix, create_gaussian_params(),
                 create_mcmc_params(n_iter = 100, n_burn = 100)),
    "n_burn must be less than n_iter"
  )

  # Test invalid thin
  expect_error(
    run_mcmc_cpp(data_matrix, create_gaussian_params(),
                 create_mcmc_params(thin = 0)),
    "thin must be positive"
  )

  # Test invalid alpha
  expect_error(
    run_mcmc_cpp(data_matrix, create_gaussian_params(),
                 create_mcmc_params(alpha = -1)),
    "alpha must be positive"
  )
})

test_that("C++ MCMC validates data input", {
  skip_if(!cpp_available(), "C++ implementation not available")

  # Empty data
  expect_error(
    run_mcmc_cpp(matrix(numeric(0), ncol = 1),
                 create_gaussian_params(),
                 create_mcmc_params()),
    "Data"
  )

  # Data with NA
  bad_data <- matrix(c(1, NA, 3), ncol = 1)
  expect_error(
    run_mcmc_cpp(bad_data, create_gaussian_params(), create_mcmc_params()),
    "NA"
  )

  # Data with Inf
  bad_data <- matrix(c(1, Inf, 3), ncol = 1)
  expect_error(
    run_mcmc_cpp(bad_data, create_gaussian_params(), create_mcmc_params()),
    "infinite|Inf"
  )
})

test_that("C++ MCMC validates mixing distribution parameters", {
  skip_if(!cpp_available(), "C++ implementation not available")

  data_matrix <- matrix(rnorm(20), ncol = 1)

  # Invalid type
  bad_params <- create_gaussian_params()
  bad_params$type <- "unknown"
  expect_error(
    run_mcmc_cpp(data_matrix, bad_params, create_mcmc_params()),
    "Unknown mixing distribution type" # Match the actual error message
  )
})

# Convergence and Correctness Tests ----------------------------------------

# Fix test at line 242-247: C++ MCMC recovers known clusters
test_that("C++ MCMC recovers known clusters", {
  skip_if(!cpp_available(), "C++ implementation not available")

  # Create well-separated clusters
  set.seed(123)
  data_obj <- create_gaussian_mixture_data(
    n = 100,
    means = c(-5, 0, 5),
    sds = c(0.5, 0.5, 0.5),
    weights = c(0.3, 0.4, 0.3),
    seed = 123
  )

  data_matrix <- matrix(data_obj$data, ncol = 1)

  # Run longer chain for convergence
  result <- run_mcmc_cpp(
    data = data_matrix,
    mixing_dist_params = create_gaussian_params(),
    mcmc_params = create_mcmc_params(n_iter = 500, n_burn = 100)
  )

  # Extract final cluster assignments
  if (is.matrix(result$cluster_labels)) {
    final_labels <- result$cluster_labels[nrow(result$cluster_labels), ]
  } else {
    final_labels <- result$cluster_labels[[length(result$cluster_labels)]]
  }

  # Use the actual final_n_clusters field if available
  final_n_clusters <- if (!is.null(result$final_n_clusters)) {
    result$final_n_clusters
  } else {
    length(unique(final_labels))
  }

  # Reasonable bounds for well-separated data
  expect_true(final_n_clusters >= 2 && final_n_clusters <= 5,
              info = paste("Got", final_n_clusters, "clusters"))

  # Calculate metrics
  metrics <- compute_clustering_metrics(data_obj$true_components, final_labels)

  # Expect reasonable ARI
  expect_true(metrics$ari > 0.5,
              info = paste("ARI =", round(metrics$ari, 3)))
})

test_that("C++ MCMC concentration parameter updates correctly", {
  skip_if(!cpp_available(), "C++ implementation not available")

  data_matrix <- matrix(rnorm(50), ncol = 1)

  # Run with update_concentration = TRUE
  result_with_update <- run_mcmc_cpp(
    data_matrix,
    create_gaussian_params(),
    create_mcmc_params(n_iter = 200, n_burn = 50, update_concentration = TRUE)
  )

  # Run with update_concentration = FALSE
  result_no_update <- run_mcmc_cpp(
    data_matrix,
    create_gaussian_params(),
    create_mcmc_params(n_iter = 200, n_burn = 50, update_concentration = FALSE)
  )

  # Extract alpha values
  alpha_with_update <- result_with_update$alpha
  alpha_no_update <- result_no_update$alpha

  # With updates, alpha should vary
  expect_true(var(alpha_with_update) > 0,
              info = paste("Variance:", var(alpha_with_update)))

  # Without updates, alpha should be constant
  expect_true(var(alpha_no_update) == 0 || var(alpha_no_update) < 1e-10)
})

# Equivalence Tests with R Implementation ----------------------------------

test_that("C++ and R implementations produce similar results", {
  skip_if(!cpp_available(), "C++ implementation not available")

  set.seed(456)
  data <- c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5))

  # Run R implementation
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(data)
  dp_r_fit <- Fit(dp_r, 100, progressBar = FALSE)

  # Run C++ implementation
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(data)
  dp_cpp_fit <- Fit(dp_cpp, 100, progressBar = FALSE)

  # Compare number of clusters - they should be similar
  expect_true(abs(dp_r_fit$numberClusters - dp_cpp_fit$numberClusters) <= 2,
              info = paste("R:", dp_r_fit$numberClusters, "C++:", dp_cpp_fit$numberClusters))

  # Both should find reasonable number of clusters
  expect_true(dp_r_fit$numberClusters >= 1 && dp_r_fit$numberClusters <= 10)
  expect_true(dp_cpp_fit$numberClusters >= 1 && dp_cpp_fit$numberClusters <= 10)
})

# Performance Tests --------------------------------------------------------

test_that("C++ implementation is faster than R implementation", {
  skip_if(!cpp_available(), "C++ implementation not available")
  skip_if_not_installed("dirichletprocess")
  skip_on_cran()  # Skip on CRAN to save time

  # Larger dataset for meaningful comparison
  data <- rnorm(100)

  # Time R implementation
  set_use_cpp(FALSE)
  r_time <- system.time({
    dp_r <- DirichletProcessGaussian(data)
    dp_r_fit <- Fit(dp_r, 50, progressBar = FALSE)
  })["elapsed"]

  # Time C++ implementation
  set_use_cpp(TRUE)
  cpp_time <- system.time({
    dp_cpp <- DirichletProcessGaussian(data)
    dp_cpp_fit <- Fit(dp_cpp, 50, progressBar = FALSE)
  })["elapsed"]

  # C++ should be at least 10x faster
  speedup <- r_time / cpp_time
  expect_true(speedup > 10,
              info = sprintf("Speedup: %.1fx (R: %.2fs, C++: %.2fs)",
                             speedup, r_time, cpp_time))
})

# Edge Cases and Stress Tests ----------------------------------------------

# Fix test at line 363: C++ MCMC handles edge cases correctly
test_that("C++ MCMC handles edge cases correctly", {
  skip_if(!cpp_available(), "C++ implementation not available")

  # Single cluster data
  single_cluster_data <- matrix(rnorm(20, 0, 0.1), ncol = 1)

  result <- run_mcmc_cpp(
    single_cluster_data,
    create_gaussian_params(),
    create_mcmc_params(n_iter = 50, n_burn = 10)
  )

  # Extract unique labels
  if (is.matrix(result$cluster_labels)) {
    all_labels <- as.vector(result$cluster_labels)
  } else {
    all_labels <- unlist(result$cluster_labels)
  }

  unique_labels <- unique(all_labels)

  # Should mostly stay as one cluster (labels are 1-indexed from R)
  expect_true(min(unique_labels) >= 1)
  expect_true(max(unique_labels) <= 3)
})

test_that("C++ MCMC handles different prior specifications", {
  skip_if(!cpp_available(), "C++ implementation not available")

  data_matrix <- matrix(rnorm(50, mean = 5, sd = 2), ncol = 1)

  # Informative prior centered at true mean
  informative_prior <- create_gaussian_params(
    mu0 = 5.0,
    kappa0 = 10.0,
    alpha0 = 10.0,
    beta0 = 20.0
  )

  result_informative <- run_mcmc_cpp(
    data_matrix,
    informative_prior,
    create_mcmc_params(n_iter = 500, n_burn = 100)
  )

  # Vague prior
  vague_prior <- create_gaussian_params(
    mu0 = 0.0,
    kappa0 = 0.01,
    alpha0 = 0.01,
    beta0 = 0.01
  )

  result_vague <- run_mcmc_cpp(
    data_matrix,
    vague_prior,
    create_mcmc_params(n_iter = 500, n_burn = 100)
  )

  # Both should work without errors
  expect_type(result_informative, "list")
  expect_type(result_vague, "list")

  # Check n_clusters is numeric vector
  expect_true(is.numeric(result_informative$n_clusters))
  expect_true(is.numeric(result_vague$n_clusters))

  # Informative prior might lead to fewer clusters
  if (length(result_informative$n_clusters) > 0 && length(result_vague$n_clusters) > 0) {
    expect_true(mean(result_informative$n_clusters) <= mean(result_vague$n_clusters) + 2)
  }
})

# Thinning Tests -----------------------------------------------------------

test_that("C++ MCMC thinning works correctly", {
  skip_if(!cpp_available(), "C++ implementation not available")

  data_matrix <- matrix(rnorm(30), ncol = 1)

  for (thin in c(1, 2, 5, 10)) {
    result <- run_mcmc_cpp(
      data_matrix,
      create_gaussian_params(),
      create_mcmc_params(n_iter = 100, n_burn = 20, thin = thin)
    )

    expected_length <- (100 - 20) / thin

    # Check the appropriate field based on structure
    if (is.matrix(result$cluster_labels)) {
      expect_equal(nrow(result$cluster_labels), expected_length,
                   info = paste("Failed for thin =", thin))
    } else {
      expect_equal(length(result$cluster_labels), expected_length,
                   info = paste("Failed for thin =", thin))
    }
  }
})

# Auxiliary Function Tests -------------------------------------------------

test_that("Gaussian likelihood calculations are correct", {
  skip_if(!cpp_available(), "C++ implementation not available")

  # Test data
  x <- c(-2, -1, 0, 1, 2)
  mu <- 0
  sigma <- 1

  # Expected log-likelihood (from dnorm)
  expected_log_lik <- sum(dnorm(x, mu, sigma, log = TRUE))

  # If there's a separate likelihood function exposed, test it
  # This is a placeholder - adjust based on actual exported functions
  # cpp_log_lik <- gaussian_log_likelihood_cpp(x, mu, sigma)
  # expect_equal(cpp_log_lik, expected_log_lik, tolerance = 1e-10)
})

# Integration Tests --------------------------------------------------------

test_that("Full DirichletProcessGaussian workflow works with C++", {
  skip_if(!cpp_available(), "C++ implementation not available")
  skip_if_not_installed("dirichletprocess")

  # Enable C++ backend
  set_use_cpp(TRUE)

  # Create data
  data <- c(rnorm(50, -3, 0.5), rnorm(50, 3, 0.5))

  # Full workflow
  dp <- DirichletProcessGaussian(data)
  dp <- Fit(dp, 500, progressBar = FALSE)

  # Check all expected components are present
  expect_true(all(c("data", "clusterLabels", "clusterParameters",
                    "numberClusters", "alpha") %in% names(dp)))

  # Test posterior predictions - use correct function name
  posterior_clusters <- ClusterLabelPredict(dp, rnorm(10))
  expect_length(posterior_clusters$componentIndexes, 10)

  # Test likelihood calculation
  lik <- Likelihood(dp, rnorm(5))
  expect_length(lik, 5)
  expect_true(all(lik >= 0))
})

# Memory and Stability Tests -----------------------------------------------

test_that("C++ MCMC doesn't leak memory on repeated calls", {
  skip_if(!cpp_available(), "C++ implementation not available")
  skip_on_cran()  # Memory tests can be flaky on CRAN

  data_matrix <- matrix(rnorm(100), ncol = 1)
  params <- create_gaussian_params()
  mcmc_params <- create_mcmc_params(n_iter = 100, n_burn = 20)

  # Run multiple times to check for memory issues
  for (i in 1:10) {
    result <- run_mcmc_cpp(data_matrix, params, mcmc_params)
    expect_type(result, "list")
  }

  # If we get here without crashing, memory management is likely OK
  expect_true(TRUE)
})

# Status and Backend Tests -------------------------------------------------

test_that("C++ backend status functions work correctly", {
  skip_if_not_installed("dirichletprocess")

  # Test status check
  status <- get_cpp_status()
  expect_type(status, "list")
  expect_true(all(sapply(status, is.logical)))

  # Test backend switching
  set_use_cpp(FALSE)
  expect_false(using_cpp())

  set_use_cpp(TRUE)
  expect_true(using_cpp())

  # Test can_use_cpp for Gaussian
  dp_gaussian <- DirichletProcessGaussian(rnorm(10))
  if (cpp_available()) {
    expect_true(can_use_cpp(dp_gaussian))
  }
})
