# tests/testthat/test_beta_distribution_complete.R
context("Beta Distribution C++ Complete Test Suite")

# Helper function to check if C++ implementation is available
skip_if_no_cpp <- function() {
  skip_if_not(exists("beta_prior_draw_cpp"), "C++ implementation not available")
}

# ===== SECTION 1: Basic Beta Distribution Tests =====

test_that("Beta prior draw C++ implementation works correctly", {
  skip_if_no_cpp()

  set.seed(123)
  priorParams <- c(2, 8)
  maxT <- 1
  n <- 100

  # Test prior draw
  result <- beta_prior_draw_cpp(priorParams, maxT, n)

  expect_equal(names(result), c("mu", "nu"))
  expect_equal(dim(result$mu), c(1, 1, n))
  expect_equal(dim(result$nu), c(1, 1, n))

  # Check all values are in valid range
  expect_true(all(result$mu >= 0 & result$mu <= maxT))
  expect_true(all(result$nu > 0))

  # Test single draw
  single_result <- beta_prior_draw_cpp(priorParams, maxT, 1)
  expect_equal(length(single_result$mu), 1)
  expect_equal(length(single_result$nu), 1)

  # Test statistical properties
  expect_true(mean(result$mu) > 0.2 && mean(result$mu) < 0.8)
  expect_true(mean(result$nu) > 0)
})

test_that("Beta likelihood C++ calculation is correct", {
  skip_if_no_cpp()

  x <- seq(0.1, 0.9, by = 0.1)
  mu <- 0.5
  nu <- 4.0
  maxT <- 1

  lik <- beta_likelihood_cpp(x, mu, nu, maxT)

  # Check dimensions
  expect_equal(length(lik), length(x))

  # Check all values are positive
  expect_true(all(lik > 0))

  # Check specific value
  a <- (mu * nu) / maxT
  b <- (1.0 - mu/maxT) * nu
  expected_lik_05 <- (1.0/maxT) * dbeta(0.5/maxT, a, b)
  expect_equal(lik[5], expected_lik_05, tolerance = 1e-10)

  # Test edge cases
  x_edge <- c(0.0, 1.0, -0.1, 1.1)
  lik_edge <- beta_likelihood_cpp(x_edge, mu, nu, maxT)
  expect_true(all(lik_edge == 1e-300))

  # Test with invalid parameters
  lik_invalid <- beta_likelihood_cpp(0.5, mu, 0.0, maxT)
  expect_equal(lik_invalid, 1e-300)
})

test_that("Beta prior density C++ calculation is correct", {
  skip_if_no_cpp()

  priorParams <- c(2, 8)
  maxT <- 1

  # Test cases
  test_cases <- list(
    list(mu = 0.5, nu = 2.0),
    list(mu = 0.2, nu = 10.0),
    list(mu = 0.8, nu = 0.5)
  )

  for (case in test_cases) {
    density <- beta_prior_density_cpp(case$mu, case$nu, priorParams, maxT)

    # Check density is positive
    expect_true(density > 0)

    # Calculate expected density
    mu_density <- 1.0 / maxT
    gamma_shape <- priorParams[1]
    gamma_rate <- priorParams[2]
    nu_density <- dgamma(1.0/case$nu, shape = gamma_shape, rate = gamma_rate) / (case$nu^2)
    expected_density <- mu_density * nu_density

    expect_equal(density, expected_density, tolerance = 1e-10)
  }

  # Test edge cases
  density_edge1 <- beta_prior_density_cpp(0.0, 2.0, priorParams, maxT)
  expect_equal(density_edge1, 1e-10)

  density_edge2 <- beta_prior_density_cpp(0.5, 1e-11, priorParams, maxT)
  expect_true(density_edge2 < 1e-5)
})

test_that("Beta posterior draw C++ produces valid samples", {
  skip_if_no_cpp()

  set.seed(456)
  priorParams <- c(2, 8)
  maxT <- 1
  mhStepSize <- c(0.1, 0.1)

  # Generate test data from known Beta
  a <- 3
  b <- 7
  x <- matrix(rbeta(30, a, b) * maxT, ncol = 1)

  # Test posterior draw
  result <- beta_posterior_draw_cpp(priorParams, maxT, mhStepSize, x,
                                    n = 5, mhDrawsVal = 250)

  expect_equal(names(result), c("mu", "nu"))
  expect_equal(dim(result$mu), c(1, 1, 5))
  expect_equal(dim(result$nu), c(1, 1, 5))

  # Check validity
  expect_true(all(result$mu > 0 & result$mu < maxT))
  expect_true(all(result$nu > 0))

  # Check that posterior mean is reasonable given the data
  true_mean <- a / (a + b) * maxT
  posterior_mu_mean <- mean(result$mu)
  expect_true(abs(posterior_mu_mean - true_mean) < 0.2)
})

# ===== SECTION 2: MCMC Implementation Tests =====

test_that("Complete Beta DP MCMC workflow works correctly", {
  skip_if_no_cpp()

  set.seed(789)

  # Generate test data from two Beta clusters
  n <- 60
  y <- c(
    rbeta(30, 2, 8),  # Cluster 1: low mean
    rbeta(30, 8, 2)   # Cluster 2: high mean
  )
  y <- sample(y)  # Shuffle

  # Create Beta DP object
  dp <- DirichletProcessBeta(y, maxT = 1, verbose = FALSE, mhDraws = 50)

  # Test initial state
  expect_equal(dp$n, n)
  expect_equal(length(dp$data), n)
  expect_equal(dp$numberClusters, 1)  # Starts with one cluster
  expect_equal(sum(dp$pointsPerCluster), n)

  # Run MCMC iterations
  dp_fitted <- Fit(dp, its = 100, progressBar = FALSE)

  # Check final state validity
  expect_true(dp_fitted$numberClusters >= 1)
  expect_true(dp_fitted$numberClusters <= n/2)  # Shouldn't have too many clusters
  expect_equal(sum(dp_fitted$pointsPerCluster), n)
  expect_true(all(dp_fitted$clusterLabels >= 1))  # R uses 1-indexing
  expect_true(all(dp_fitted$clusterLabels <= dp_fitted$numberClusters))

  # Check that parameters are reasonable
  mu_params <- dp_fitted$clusterParameters$mu
  nu_params <- dp_fitted$clusterParameters$nu

  expect_true(all(mu_params > 0 & mu_params < 1))
  expect_true(all(nu_params > 0))

  # Should find approximately 2 clusters
  expect_true(dp_fitted$numberClusters >= 1 && dp_fitted$numberClusters <= 4)

  # Check alpha chain
  expect_true(length(dp_fitted$alphaChain) == 100)
  expect_true(all(dp_fitted$alphaChain > 0))
})

test_that("Beta DP handles different data patterns correctly", {
  skip_if_no_cpp()

  # Test 1: Single cluster data
  set.seed(111)
  y_single <- rbeta(50, 5, 5)  # Symmetric, single mode
  dp_single <- DirichletProcessBeta(y_single, verbose = FALSE)
  dp_single <- Fit(dp_single, its = 50, progressBar = FALSE)

  expect_true(dp_single$numberClusters >= 1 && dp_single$numberClusters <= 3)

  # Test 2: Well-separated clusters
  set.seed(222)
  y_separated <- c(
    rbeta(25, 1, 10),  # Very low values
    rbeta(25, 10, 1)   # Very high values
  )
  dp_separated <- DirichletProcessBeta(y_separated, verbose = FALSE)
  dp_separated <- Fit(dp_separated, its = 50, progressBar = FALSE)

  # Should find 2 clusters with high probability
  expect_true(dp_separated$numberClusters >= 2)

  # Test 3: Small dataset
  set.seed(333)
  y_small <- rbeta(10, 3, 3)
  dp_small <- DirichletProcessBeta(y_small, verbose = FALSE)
  dp_small <- Fit(dp_small, its = 50, progressBar = FALSE)

  expect_true(dp_small$numberClusters >= 1 && dp_small$numberClusters <= 5)
  expect_equal(sum(dp_small$pointsPerCluster), 10)
})

test_that("Beta DP parameter updates are sensible", {
  skip_if_no_cpp()

  set.seed(444)

  # Generate data with known parameters
  true_mu1 <- 0.2
  true_mu2 <- 0.8
  true_nu <- 10

  y <- c(
    rbeta(40, true_mu1 * true_nu, (1 - true_mu1) * true_nu),
    rbeta(40, true_mu2 * true_nu, (1 - true_mu2) * true_nu)
  )

  dp <- DirichletProcessBeta(y, verbose = FALSE, mhStepSize = c(0.05, 0.05))
  dp <- Fit(dp, its = 200, progressBar = FALSE)

  # Check that we found clusters near the true values
  mu_params <- sort(dp$clusterParameters$mu)

  if (dp$numberClusters >= 2) {
    # Find the two main clusters (with most points)
    points_per_cluster <- dp$pointsPerCluster
    top_clusters <- order(points_per_cluster, decreasing = TRUE)[1:2]
    main_mus <- sort(dp$clusterParameters$mu[top_clusters])

    # Check they're close to true values
    expect_true(abs(main_mus[1] - true_mu1) < 0.1)
    expect_true(abs(main_mus[2] - true_mu2) < 0.1)
  }
})

test_that("Beta DP alpha updates follow expected behavior", {
  skip_if_no_cpp()

  set.seed(555)

  # Test with different alpha priors
  y <- rbeta(50, 3, 3)

  # High alpha prior - expect more clusters
  dp_high <- DirichletProcessBeta(y, alphaPriors = c(10, 2), verbose = FALSE)
  dp_high <- Fit(dp_high, its = 100, progressBar = FALSE)

  # Low alpha prior - expect fewer clusters
  dp_low <- DirichletProcessBeta(y, alphaPriors = c(1, 10), verbose = FALSE)
  dp_low <- Fit(dp_low, its = 100, progressBar = FALSE)

  # Average number of clusters should reflect alpha
  mean_alpha_high <- mean(tail(dp_high$alphaChain, 50))
  mean_alpha_low <- mean(tail(dp_low$alphaChain, 50))

  expect_true(mean_alpha_high > mean_alpha_low)
})

test_that("Beta DP handles edge cases gracefully", {
  skip_if_no_cpp()

  # Test 1: Extreme data values
  set.seed(666)
  y_extreme <- c(
    rep(0.001, 5),  # Very small
    rep(0.999, 5)   # Very large
  )

  expect_error({
    dp_extreme <- DirichletProcessBeta(y_extreme, verbose = FALSE)
    dp_extreme <- Fit(dp_extreme, its = 20, progressBar = FALSE)
  }, NA)  # Should not error

  # Test 2: Single data point
  y_single <- 0.5
  dp_single <- DirichletProcessBeta(y_single, verbose = FALSE)
  dp_single <- Fit(dp_single, its = 10, progressBar = FALSE)

  expect_equal(dp_single$numberClusters, 1)
  expect_equal(dp_single$pointsPerCluster, 1)

  # Test 3: Identical data points
  y_identical <- rep(0.7, 20)
  dp_identical <- DirichletProcessBeta(y_identical, verbose = FALSE)
  dp_identical <- Fit(dp_identical, its = 30, progressBar = FALSE)

  # Should find very few clusters
  expect_true(dp_identical$numberClusters <= 3)
})

# ===== SECTION 3: Integration and Convergence Tests =====

test_that("Beta DP MCMC chain shows convergence", {
  skip_if_no_cpp()

  set.seed(777)

  # Generate data
  y <- c(rbeta(30, 2, 8), rbeta(30, 8, 2))
  dp <- DirichletProcessBeta(y, verbose = FALSE)

  # Run longer chain
  dp <- Fit(dp, its = 500, progressBar = FALSE)

  # Check convergence of number of clusters
  cluster_chain <- dp$weightsChain
  n_clusters <- apply(cluster_chain > 0, 2, sum)

  # Compare first and last halves
  first_half <- n_clusters[1:250]
  second_half <- n_clusters[251:500]

  # Means should be similar if converged
  expect_true(abs(mean(first_half) - mean(second_half)) < 0.5)

  # Check alpha convergence
  alpha_first <- dp$alphaChain[1:250]
  alpha_second <- dp$alphaChain[251:500]

  expect_true(abs(mean(alpha_first) - mean(alpha_second)) <
                0.3 * mean(dp$alphaChain))

  # Variance should stabilize
  expect_true(var(alpha_second) < 2 * var(alpha_first))
})

test_that("Beta DP produces consistent results with same seed", {
  skip_if_no_cpp()

  y <- rbeta(40, 3, 7)

  # Run 1
  set.seed(888)
  dp1 <- DirichletProcessBeta(y, verbose = FALSE)
  dp1 <- Fit(dp1, its = 50, progressBar = FALSE)

  # Run 2 with same seed
  set.seed(888)
  dp2 <- DirichletProcessBeta(y, verbose = FALSE)
  dp2 <- Fit(dp2, its = 50, progressBar = FALSE)

  # Results should be identical
  expect_equal(dp1$numberClusters, dp2$numberClusters)
  expect_equal(dp1$clusterLabels, dp2$clusterLabels)
  expect_equal(dp1$alpha, dp2$alpha)
  expect_equal(dp1$clusterParameters$mu, dp2$clusterParameters$mu)
  expect_equal(dp1$clusterParameters$nu, dp2$clusterParameters$nu)
})

test_that("Beta DP likelihood calculations are correct", {
  skip_if_no_cpp()

  set.seed(999)

  y <- rbeta(30, 4, 6)
  dp <- DirichletProcessBeta(y, verbose = FALSE)
  dp <- Fit(dp, its = 100, progressBar = FALSE)

  # Calculate likelihood for the final state
  final_lik <- LikelihoodDP(dp)

  expect_true(is.numeric(final_lik))
  expect_true(is.finite(final_lik))
  expect_true(final_lik < 0)  # Log likelihood should be negative

  # Check likelihood chain
  expect_equal(length(dp$likelihoodChain), 100)
  expect_true(all(is.finite(dp$likelihoodChain)))

  # Likelihood should generally increase (with some randomness)
  smooth_lik <- filter(dp$likelihoodChain, rep(1/10, 10), sides = 1)
  smooth_lik <- smooth_lik[!is.na(smooth_lik)]
  expect_true(tail(smooth_lik, 1) > head(smooth_lik, 1))
})

test_that("Beta DP posterior sampling works", {
  skip_if_no_cpp()

  set.seed(1234)

  y <- rbeta(40, 3, 3)
  dp <- DirichletProcessBeta(y, verbose = FALSE)
  dp <- Fit(dp, its = 100, progressBar = FALSE)

  # Draw from posterior
  posterior_sample <- PosteriorFunction(dp, 100)

  expect_true(is.function(posterior_sample))

  # Test posterior function
  test_points <- seq(0.1, 0.9, by = 0.1)
  posterior_values <- posterior_sample(test_points)

  expect_equal(length(posterior_values), length(test_points))
  expect_true(all(posterior_values >= 0))
  expect_true(all(is.finite(posterior_values)))

  # Posterior should integrate to approximately 1
  integrate_result <- integrate(posterior_sample, 0, 1)
  expect_true(abs(integrate_result$value - 1) < 0.1)
})

test_that("Beta DP cluster assignment predictions work", {
  skip_if_no_cpp()

  set.seed(5678)

  # Train on subset
  y_train <- c(rbeta(20, 2, 8), rbeta(20, 8, 2))
  dp <- DirichletProcessBeta(y_train, verbose = FALSE)
  dp <- Fit(dp, its = 100, progressBar = FALSE)

  # Predict on new data
  y_test <- c(0.1, 0.2, 0.8, 0.9)  # Should assign to different clusters
  pred_clusters <- ClusterLabelPredict(dp, y_test)

  expect_equal(length(pred_clusters), length(y_test))
  expect_true(all(pred_clusters >= 1))
  expect_true(all(pred_clusters <= dp$numberClusters + 1))  # Can create new cluster

  # Low values should be in same cluster, high values in same cluster
  expect_equal(pred_clusters[1], pred_clusters[2])
  expect_equal(pred_clusters[3], pred_clusters[4])
  expect_true(pred_clusters[1] != pred_clusters[3])  # Different clusters
})

test_that("Beta DP methods handle various maxT values", {
  skip_if_no_cpp()

  set.seed(9999)

  # Test with different maxT values
  maxT_values <- c(1, 10, 100)

  for (maxT in maxT_values) {
    y <- rbeta(30, 3, 7) * maxT

    dp <- DirichletProcessBeta(y, maxT = maxT, verbose = FALSE)
    dp <- Fit(dp, its = 50, progressBar = FALSE)

    # Check parameters are scaled correctly
    expect_true(all(dp$clusterParameters$mu > 0))
    expect_true(all(dp$clusterParameters$mu < maxT))
    expect_true(all(dp$data >= 0))
    expect_true(all(dp$data <= maxT))

    # Likelihood should be finite
    lik <- LikelihoodDP(dp)
    expect_true(is.finite(lik))
  }
})

# Performance test (only run if explicitly requested)
test_that("Beta DP performance is reasonable", {
  skip_on_cran()
  skip_if_not(interactive(), "Performance test only run interactively")

  set.seed(1111)

  # Test scaling with data size
  n_values <- c(50, 100, 200)
  times <- numeric(length(n_values))

  for (i in seq_along(n_values)) {
    n <- n_values[i]
    y <- rbeta(n, 3, 7)

    time_start <- Sys.time()
    dp <- DirichletProcessBeta(y, verbose = FALSE)
    dp <- Fit(dp, its = 100, progressBar = FALSE)
    time_end <- Sys.time()

    times[i] <- as.numeric(time_end - time_start, units = "secs")
  }

  cat("\nBeta DP Performance (100 iterations):\n")
  for (i in seq_along(n_values)) {
    cat(sprintf("n = %d: %.3f seconds\n", n_values[i], times[i]))
  }

  # Time should scale roughly linearly with n
  expect_true(times[3] < times[1] * 5)  # Not more than 5x slower for 4x data
})
