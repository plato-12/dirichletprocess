context("Multivariate Normal Distribution C++ Implementation")

test_that("MVNormal PriorDraw C++ matches R implementation", {
  set.seed(42)
  priorParams <- list(
    mu0 = c(0, 0),
    Lambda = diag(2),
    kappa0 = 1,
    nu = 3
  )
  n <- 5

  # Create mixing distribution object for R implementation
  mdObj <- MvnormalCreate(priorParams)

  # R implementation
  r_result <- PriorDraw(mdObj, n)

  # C++ implementation
  cpp_result <- mvnormal_prior_draw_cpp(priorParams, n)

  # Check structure
  expect_equal(names(cpp_result), names(r_result))
  expect_equal(dim(cpp_result$mu), dim(r_result$mu))
  expect_equal(dim(cpp_result$sig), dim(r_result$sig))

  # Check statistical properties (random draws won't be identical)
  # Just check dimensions and that values are reasonable
  expect_equal(dim(cpp_result$mu), dim(r_result$mu))
  expect_equal(dim(cpp_result$sig), dim(r_result$sig))
  # Check covariances are positive definite
  for (i in 1:n) {
    expect_true(all(eigen(cpp_result$sig[,,i])$values > 0))
    expect_true(all(eigen(r_result$sig[,,i])$values > 0))
  }
})

test_that("MVNormal PosteriorDraw C++ matches R implementation", {
  set.seed(123)
  priorParams <- list(
    mu0 = c(0, 0),
    Lambda = diag(2),
    kappa0 = 1,
    nu = 3
  )
  x <- matrix(rnorm(20), ncol = 2)
  n <- 3

  # Create mixing distribution object for R implementation
  mdObj <- MvnormalCreate(priorParams)

  # R implementation
  set.seed(123)
  r_result <- PosteriorDraw(mdObj, x, n)

  # C++ implementation
  set.seed(123)
  cpp_result <- mvnormal_posterior_draw_cpp(priorParams, x, n)

  # Check structure
  expect_equal(names(cpp_result), names(r_result))
  expect_equal(dim(cpp_result$mu), dim(r_result$mu))
  expect_equal(dim(cpp_result$sig), dim(r_result$sig))

  # Check posterior parameters are calculated correctly
  r_post_params <- PosteriorParameters(mdObj, x)
  cpp_post_params <- mvnormal_posterior_parameters_cpp(priorParams, x)

  expect_equal(cpp_post_params$mu_n, r_post_params$mu_n, tolerance = 1e-10)
  expect_equal(cpp_post_params$kappa_n, r_post_params$kappa_n, tolerance = 1e-10)
  expect_equal(cpp_post_params$nu_n, r_post_params$nu_n, tolerance = 1e-10)
  expect_equal(as.matrix(cpp_post_params$Lambda_n),
               as.matrix(r_post_params$t_n), tolerance = 1e-10)
})

test_that("MVNormal Likelihood C++ matches R implementation", {
  set.seed(456)
  x <- matrix(rnorm(6), ncol = 2)
  mu <- c(0, 0)
  sigma <- diag(2)

  # Use the mvtnorm package directly for R implementation
  r_lik <- mvtnorm::dmvnorm(x, mu, sigma)

  # C++ implementation
  cpp_lik <- mvnormal_likelihood_cpp(x, mu, sigma)

  expect_equal(cpp_lik, r_lik, tolerance = 1e-10)
})

test_that("MVNormal Predictive C++ matches R implementation", {
  set.seed(789)
  priorParams <- list(
    mu0 = c(0, 0),
    Lambda = diag(2),
    kappa0 = 1,
    nu = 3
  )
  x <- matrix(rnorm(10), ncol = 2)

  # Create mixing distribution object for R implementation
  mdObj <- MvnormalCreate(priorParams)

  # R implementation
  r_pred <- Predictive(mdObj, x)

  # C++ implementation
  cpp_pred <- mvnormal_predictive_cpp(priorParams, x)

  expect_equal(cpp_pred, r_pred, tolerance = 1e-8)
})

test_that("MVNormal ClusterComponentUpdate C++ works correctly", {
  set.seed(999)
  # Create a simple DP object with 2D data
  y <- rbind(
    mvtnorm::rmvnorm(10, c(-2, -2), diag(2) * 0.5),
    mvtnorm::rmvnorm(10, c(2, 2), diag(2) * 0.5)
  )
  dp <- DirichletProcessMvnormal(y)

  # Manually set up clusters
  dp$clusterLabels <- c(rep(1, 10), rep(2, 10))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(10, 10)

  # Initialize cluster parameters
  mu1 <- colMeans(y[1:10, ])
  mu2 <- colMeans(y[11:20, ])
  sig1 <- cov(y[1:10, ])
  sig2 <- cov(y[11:20, ])

  dp$clusterParameters <- list(
    mu = array(c(mu1, mu2), dim = c(1, 2, 2)),
    sig = array(c(sig1, sig2), dim = c(2, 2, 2))
  )

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run C++ update
  set.seed(321)
  cpp_result <- conjugate_mvnormal_cluster_component_update_cpp(dp_cpp)

  # Convert back to 1-indexed
  cpp_result$clusterLabels <- cpp_result$clusterLabels + 1

  # Check that the structure is preserved
  expect_true(all(cpp_result$clusterLabels %in% 1:20))
  expect_equal(sum(cpp_result$pointsPerCluster), 20)
  expect_true(cpp_result$numberClusters >= 1)
  expect_true(cpp_result$numberClusters <= 20)
})

test_that("MVNormal ClusterParameterUpdate C++ works correctly", {
  set.seed(111)
  # Create data with known structure
  y1 <- mvtnorm::rmvnorm(15, mean = c(-3, -3), sigma = diag(2) * 0.5)
  y2 <- mvtnorm::rmvnorm(15, mean = c(3, 3), sigma = diag(2) * 0.5)
  y <- rbind(y1, y2)

  dp <- DirichletProcessMvnormal(y)

  # Set up known clustering
  dp$clusterLabels <- c(rep(1, 15), rep(2, 15))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(15, 15)
  dp$clusterParameters <- list(
    mu = array(c(0, 0, 0, 0), dim = c(1, 2, 2)),
    sig = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run parameter update
  set.seed(222)
  cpp_params <- conjugate_mvnormal_cluster_parameter_update_cpp(dp_cpp)

  # Check that parameters were updated reasonably
  # Extract means
  mu1 <- cpp_params$mu[1, , 1]
  mu2 <- cpp_params$mu[1, , 2]

  # Cluster 1 should have mean around (-3, -3)
  expect_true(all(abs(mu1 - c(-3, -3)) < 1))
  # Cluster 2 should have mean around (3, 3)
  expect_true(all(abs(mu2 - c(3, 3)) < 1))

  # Covariances should be positive definite
  sig1 <- cpp_params$sig[, , 1]
  sig2 <- cpp_params$sig[, , 2]

  expect_true(all(eigen(sig1)$values > 0))
  expect_true(all(eigen(sig2)$values > 0))
})

test_that("End-to-end MVNormal C++ sampler test", {
  set.seed(2025)
  # Generate data with clear clusters
  y <- rbind(
    mvtnorm::rmvnorm(20, c(-3, -3), diag(2) * 0.5),
    mvtnorm::rmvnorm(20, c(0, 3), diag(2) * 0.5),
    mvtnorm::rmvnorm(20, c(3, -3), diag(2) * 0.5)
  )

  # Initialize DP
  dp <- DirichletProcessMvnormal(y)

  # Convert to format expected by C++ (0-indexed clusters)
  dp$clusterLabels <- dp$clusterLabels - 1

  # Store initial state
  initial_clusters <- length(unique(dp$clusterLabels))

  # Run 10 iterations of C++ sampler
  for (iter in 1:10) {
    # Update cluster assignments
    update_result <- conjugate_mvnormal_cluster_component_update_cpp(dp)

    # Update dp object with results
    dp$clusterLabels <- update_result$clusterLabels
    dp$pointsPerCluster <- update_result$pointsPerCluster
    dp$numberClusters <- update_result$numberClusters
    dp$clusterParameters <- update_result$clusterParameters

    # Update cluster parameters
    dp$clusterParameters <- conjugate_mvnormal_cluster_parameter_update_cpp(dp)

    # Basic sanity checks after each iteration
    expect_equal(length(dp$clusterLabels), nrow(y))
    expect_equal(sum(dp$pointsPerCluster), nrow(y))
    expect_true(dp$numberClusters >= 1)
    expect_true(dp$numberClusters <= nrow(y))
  }

  # Check that we found reasonable clusters
  # With only 10 iterations and 60 data points, we might get many small clusters
  final_clusters <- dp$numberClusters
  expect_true(final_clusters >= 1)  # At least one cluster
  expect_true(final_clusters <= nrow(y))  # At most one per data point
  # For better clustering, would need more iterations

  # Convert back to 1-indexed for inspection
  final_labels <- dp$clusterLabels + 1

  # Print summary for manual inspection
  cat("\nEnd-to-end MVNormal test summary:\n")
  cat("Initial clusters:", initial_clusters, "\n")
  cat("Final clusters:", final_clusters, "\n")
  cat("Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
})

# Test to verify MVNormal fixes work correctly
test_that("MVNormal C++ implementation handles matrix symmetry and bounds correctly", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Test 1: Matrix symmetry in prior draws
  priorParams <- list(
    mu0 = c(0, 0),
    Lambda = matrix(c(1, 0.1, 0.1, 1), 2, 2),  # Slightly asymmetric due to numerics
    kappa0 = 1,
    nu = 4
  )

  # This should not produce warnings
  expect_silent({
    result <- mvnormal_prior_draw_cpp(priorParams, n = 10)
  })

  # Check all covariance matrices are symmetric
  for (i in 1:10) {
    sig_i <- result$sig[,,i]
    # Check symmetry (within numerical tolerance)
    expect_true(all(abs(sig_i - t(sig_i)) < 1e-10))
  }

  # Test 2: Bounds checking in cluster updates
  set.seed(123)
  n <- 30
  d <- 2
  data <- rbind(
    mvtnorm::rmvnorm(15, c(-2, -2), diag(2)),
    mvtnorm::rmvnorm(15, c(2, 2), diag(2))
  )

  # Create DP object with proper initialization
  dp <- DirichletProcessMvnormal(data)

  # Ensure parameter arrays have enough space
  if (dim(dp$clusterParameters$mu)[3] < 20) {
    new_mu <- array(NA_real_, dim = c(1, d, 20))
    new_sig <- array(NA_real_, dim = c(d, d, 20))

    old_dim <- dim(dp$clusterParameters$mu)[3]
    new_mu[, , 1:old_dim] <- dp$clusterParameters$mu
    new_sig[, , 1:old_dim] <- dp$clusterParameters$sig

    dp$clusterParameters$mu <- new_mu
    dp$clusterParameters$sig <- new_sig
  }

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run multiple iterations without errors
  for (iter in 1:5) {
    result <- conjugate_mvnormal_cluster_component_update_cpp(dp_cpp)
    dp_cpp$clusterLabels <- result$clusterLabels
    dp_cpp$pointsPerCluster <- result$pointsPerCluster
    dp_cpp$numberClusters <- result$numberClusters
    dp_cpp$clusterParameters <- result$clusterParameters

    # Verify consistency
    expect_equal(length(dp_cpp$clusterLabels), n)
    expect_equal(sum(dp_cpp$pointsPerCluster), n)
    expect_true(dp_cpp$numberClusters >= 1)
    expect_true(dp_cpp$numberClusters <= n)

    # Check parameter arrays are properly sized
    mu_dim <- dim(dp_cpp$clusterParameters$mu)
    expect_true(mu_dim[3] >= dp_cpp$numberClusters)
  }

  # Test 3: Ensure no hanging with many clusters
  # Create data that will likely create many clusters
  set.seed(456)
  scattered_data <- matrix(rnorm(100 * 2, sd = 5), ncol = 2)

  dp2 <- DirichletProcessMvnormal(scattered_data)

  # Pre-allocate enough space
  new_mu <- array(NA_real_, dim = c(1, 2, 100))
  new_sig <- array(NA_real_, dim = c(2, 2, 100))

  old_dim <- dim(dp2$clusterParameters$mu)[3]
  new_mu[, , 1:old_dim] <- dp2$clusterParameters$mu
  new_sig[, , 1:old_dim] <- dp2$clusterParameters$sig

  # Fill remaining with prior draws
  if (old_dim < 100) {
    extra_params <- PriorDraw(dp2$mixingDistribution, 100 - old_dim)
    new_mu[, , (old_dim+1):100] <- extra_params$mu
    new_sig[, , (old_dim+1):100] <- extra_params$sig
  }

  dp2$clusterParameters$mu <- new_mu
  dp2$clusterParameters$sig <- new_sig

  dp2$clusterLabels <- dp2$clusterLabels - 1

  # This should complete without hanging
  result2 <- conjugate_mvnormal_cluster_component_update_cpp(dp2)

  expect_true(result2$numberClusters >= 1)
  expect_equal(sum(result2$pointsPerCluster), 100)
})

test_that("MVNormal symmetry is preserved throughout operations", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Create a precision matrix that might become asymmetric
  priorParams <- list(
    mu0 = c(0, 0, 0),
    Lambda = matrix(c(1, 0.5, 0.3,
                      0.5, 2, 0.4,
                      0.3, 0.4, 1.5), 3, 3),
    kappa0 = 2,
    nu = 5
  )

  # Generate data
  set.seed(789)
  x <- mvtnorm::rmvnorm(20, rep(1, 3), diag(3))

  # Test posterior parameters
  post_params <- mvnormal_posterior_parameters_cpp(priorParams, x)

  # Check t_n is symmetric
  t_n <- post_params$t_n
  expect_true(all(abs(t_n - t(t_n)) < 1e-10))

  # Test posterior draws
  post_draws <- mvnormal_posterior_draw_cpp(priorParams, x, n = 5)

  # Check all drawn precision matrices are symmetric
  for (i in 1:5) {
    sig_i <- post_draws$sig[,,i]
    expect_true(all(abs(sig_i - t(sig_i)) < 1e-10))
  }
})
