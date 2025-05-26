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
  set.seed(42)
  r_result <- PriorDraw(mdObj, n)

  # C++ implementation
  set.seed(42)
  cpp_result <- mvnormal_prior_draw_cpp(priorParams, n)

  # Check structure
  expect_equal(names(cpp_result), names(r_result))
  expect_equal(dim(cpp_result$mu), dim(r_result$mu))
  expect_equal(dim(cpp_result$sig), dim(r_result$sig))

  # Check statistical properties
  expect_equal(mean(cpp_result$mu), mean(r_result$mu), tolerance = 0.3)
  expect_equal(mean(cpp_result$sig), mean(r_result$sig), tolerance = 0.3)
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
  final_clusters <- dp$numberClusters
  expect_true(final_clusters >= 2 && final_clusters <= 6)

  # Convert back to 1-indexed for inspection
  final_labels <- dp$clusterLabels + 1

  # Print summary for manual inspection
  cat("\nEnd-to-end MVNormal test summary:\n")
  cat("Initial clusters:", initial_clusters, "\n")
  cat("Final clusters:", final_clusters, "\n")
  cat("Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
})
