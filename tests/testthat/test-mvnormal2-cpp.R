# tests/testthat/test-mvnormal2-cpp.R

context("MVNormal2 C++ Implementation")

test_that("MVNormal2 prior draw works correctly", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Set up prior parameters
  priorParams <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2),
    phi0 = diag(2),
    nu0 = 4
  )

  # Test single draw
  set.seed(123)
  result <- mvnormal2_prior_draw_cpp(priorParams, n = 1)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_equal(names(result), c("mu", "sig"))

  # Check dimensions
  expect_equal(dim(result$mu), c(1, 2, 1))
  expect_equal(dim(result$sig), c(2, 2, 1))

  # Test multiple draws
  set.seed(123)
  result_multi <- mvnormal2_prior_draw_cpp(priorParams, n = 10)

  expect_equal(dim(result_multi$mu), c(1, 2, 10))
  expect_equal(dim(result_multi$sig), c(2, 2, 10))

  # Check that covariance matrices are positive definite
  for (i in 1:10) {
    sig_i <- result_multi$sig[,,i]
    eigenvals <- eigen(sig_i)$values
    expect_true(all(eigenvals > 0))
  }
})

test_that("MVNormal2 posterior draw works correctly", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Set up prior parameters
  priorParams <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2),
    phi0 = diag(2),
    nu0 = 4
  )

  # Generate some data
  set.seed(456)
  n_obs <- 20
  true_mu <- c(2, -1)
  true_sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  x <- mvtnorm::rmvnorm(n_obs, true_mu, true_sigma)

  # Test posterior draw
  set.seed(789)
  result <- mvnormal2_posterior_draw_cpp(priorParams, x, n = 100)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_equal(dim(result$mu), c(1, 2, 100))
  expect_equal(dim(result$sig), c(2, 2, 100))

  # Check that posterior means are reasonable
  posterior_mu_mean <- rowMeans(result$mu[1,,])
  expect_true(all(abs(posterior_mu_mean - true_mu) < 1))

  # Test with single observation
  x_single <- matrix(c(1, 2), nrow = 1)
  result_single <- mvnormal2_posterior_draw_cpp(priorParams, x_single, n = 10)
  expect_equal(dim(result_single$mu), c(1, 2, 10))

  # Test with empty data (should return prior)
  x_empty <- matrix(numeric(0), ncol = 2)
  result_empty <- mvnormal2_posterior_draw_cpp(priorParams, x_empty, n = 5)
  expect_equal(dim(result_empty$mu), c(1, 2, 5))
})

test_that("MVNormal2 likelihood calculation works correctly", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Set up test data
  x <- c(1, 2)
  theta <- list(
    mu = array(c(0, 0, 1, 1, 2, 2), dim = c(1, 2, 3)),
    sig = array(c(rep(c(1, 0, 0, 1), 3)), dim = c(2, 2, 3))
  )

  # Calculate likelihood
  result <- mvnormal2_likelihood_cpp(x, theta)

  expect_is(result, "numeric")
  expect_equal(length(result), 3)
  expect_true(all(result > 0))
  expect_true(all(result < 1))

  # Compare with R implementation
  r_result <- numeric(3)
  for (i in 1:3) {
    r_result[i] <- mvtnorm::dmvnorm(x, theta$mu[1,,i], theta$sig[,,i])
  }

  expect_equal(result, r_result, tolerance = 1e-6)
})

test_that("NonConjugate MVNormal2 cluster updates work", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Create a simple DP object
  set.seed(123)
  n <- 30
  d <- 2
  data <- rbind(
    mvtnorm::rmvnorm(15, c(-2, -2), diag(2)),
    mvtnorm::rmvnorm(15, c(2, 2), diag(2))
  )

  # Initialize DP object
  dpObj <- list(
    data = data,
    n = n,
    alpha = 1.0,
    alphaPriorParameters = c(2, 4),
    clusterLabels = rep(1:2, each = 15) - 1,  # 0-indexed
    pointsPerCluster = c(15, 15),
    numberClusters = 2,
    clusterParameters = list(
      mu = array(c(-2, -2, 2, 2), dim = c(1, 2, 2)),
      sig = array(rep(c(1, 0, 0, 1), 2), dim = c(2, 2, 2))
    ),
    mixingDistribution = list(
      priorParameters = list(
        mu0 = c(0, 0),
        sigma0 = diag(2),
        phi0 = diag(2),
        nu0 = 4
      )
    ),
    m = 3,
    mhDraws = 100
  )

  # Test cluster component update
  result <- nonconjugate_mvnormal2_cluster_component_update_cpp(dpObj)

  expect_is(result, "list")
  expect_equal(length(result$clusterLabels), n)
  expect_true(all(result$clusterLabels >= 0))
  expect_equal(sum(result$pointsPerCluster), n)
  expect_equal(length(result$pointsPerCluster), result$numberClusters)

  # Test cluster parameter update
  result_params <- nonconjugate_mvnormal2_cluster_parameter_update_cpp(dpObj)

  expect_is(result_params, "list")
  expect_equal(length(result_params), 2)
  expect_equal(names(result_params), c("mu", "sig"))
})

test_that("Hierarchical MVNormal2 creation works", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))
  skip_if_not(requireNamespace("gtools", quietly = TRUE))

  n <- 3
  priorParams <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2),
    phi0 = diag(2),
    nu0 = 4
  )
  alphaPrior <- c(2, 4)
  gammaPrior <- c(2, 2)
  num_sticks <- 10

  # Test creation
  result <- hierarchical_mvnormal2_mixing_create_cpp(
    n, priorParams, alphaPrior, gammaPrior, num_sticks
  )

  expect_is(result, "list")
  expect_equal(length(result), n)

  # Check each mixing distribution
  for (i in 1:n) {
    md <- result[[i]]
    expect_equal(md$distribution, "mvnormal2")
    expect_false(md$conjugate)
    expect_is(md$theta_k, "list")
    expect_is(md$beta_k, "numeric")
    expect_is(md$gamma, "numeric")
    expect_is(md$alpha, "numeric")
    expect_is(md$pi_k, "numeric")

    # Check dimensions
    expect_equal(length(md$beta_k), num_sticks)
    expect_equal(length(md$pi_k), num_sticks)

    # Check that weights sum to approximately 1
    expect_equal(sum(md$beta_k), 1, tolerance = 1e-10)
    expect_equal(sum(md$pi_k), 1, tolerance = 1e-10)
  }
})

test_that("Hierarchical MVNormal2 fitting works", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))
  skip_if_not(requireNamespace("gtools", quietly = TRUE))

  # Create simple hierarchical data
  set.seed(456)
  dataList <- list(
    mvtnorm::rmvnorm(20, c(-1, -1), diag(2)),
    mvtnorm::rmvnorm(25, c(1, 1), diag(2)),
    mvtnorm::rmvnorm(30, c(0, 2), diag(2))
  )

  # Create hierarchical DP
  priorParams <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2) * 2,
    phi0 = diag(2) * 2,
    nu0 = 5
  )

  hdp <- DirichletProcessHierarchicalMvnormal2(
    dataList = dataList,
    g0Priors = priorParams,
    gammaPriors = c(2, 2),
    alphaPriors = c(2, 4),
    numSticks = 20,
    numInitialClusters = 1
  )

  # Convert to format for C++
  hdp_list <- list(
    indDP = hdp$indDP,
    globalParameters = hdp$globalParameters,
    globalStick = hdp$globalStick,
    gamma = hdp$gamma,
    gammaPriors = hdp$gammaPriors
  )

  # Fit for a few iterations
  set.seed(789)
  result <- hierarchical_mvnormal2_fit_cpp(hdp_list, iterations = 5,
                                           updatePrior = FALSE,
                                           progressBar = FALSE)

  expect_is(result, "list")
  expect_equal(length(result$indDP), 3)
  expect_is(result$globalParameters, "list")
  expect_is(result$gamma, "numeric")

  # Check that each individual DP has been updated
  for (i in 1:3) {
    dp <- result$indDP[[i]]
    expect_equal(nrow(dp$data), nrow(dataList[[i]]))
    expect_true(all(dp$clusterLabels >= 0))
    expect_equal(sum(dp$pointsPerCluster), nrow(dataList[[i]]))
  }
})

test_that("MVNormal2 R and C++ implementations are consistent", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Set up test case
  set.seed(999)
  priorParams <- list(
    mu0 = c(1, -1),
    sigma0 = matrix(c(2, 0.5, 0.5, 2), 2, 2),
    phi0 = matrix(c(3, 1, 1, 3), 2, 2),
    nu0 = 6
  )

  # Test data
  x <- matrix(rnorm(20), ncol = 2)

  # Create mixing distribution objects
  mdObj_r <- Mvnormal2Create(priorParams)

  # Test prior draws
  set.seed(111)
  r_prior <- PriorDraw.mvnormal2(mdObj_r, n = 5)

  set.seed(111)
  cpp_prior <- mvnormal2_prior_draw_cpp(priorParams, n = 5)

  # Check dimensions match
  expect_equal(dim(r_prior$mu), dim(cpp_prior$mu))
  expect_equal(dim(r_prior$sig), dim(cpp_prior$sig))

  # Test posterior draws (these won't be identical due to different algorithms)
  set.seed(222)
  r_post <- PosteriorDraw.mvnormal2(mdObj_r, x, n = 10)

  set.seed(222)
  cpp_post <- mvnormal2_posterior_draw_cpp(priorParams, x, n = 10)

  # Check dimensions
  expect_equal(dim(r_post$mu), dim(cpp_post$mu))
  expect_equal(dim(r_post$sig), dim(cpp_post$sig))

  # Check that means are in reasonable range
  r_mu_mean <- apply(r_post$mu, 2, mean)
  cpp_mu_mean <- apply(cpp_post$mu, 2, mean)
  data_mean <- colMeans(x)

  expect_true(all(abs(r_mu_mean - data_mean) < 2))
  expect_true(all(abs(cpp_mu_mean - data_mean) < 2))
})

test_that("MVNormal2 handles edge cases correctly", {
  priorParams <- list(
    mu0 = c(0, 0, 0),
    sigma0 = diag(3),
    phi0 = diag(3),
    nu0 = 5
  )

  # Test with different dimensions
  result_3d <- mvnormal2_prior_draw_cpp(priorParams, n = 5)
  expect_equal(dim(result_3d$mu), c(1, 3, 5))
  expect_equal(dim(result_3d$sig), c(3, 3, 5))

  # Test with 1D (should still work)
  priorParams_1d <- list(
    mu0 = 0,
    sigma0 = matrix(1, 1, 1),
    phi0 = matrix(1, 1, 1),
    nu0 = 3
  )

  result_1d <- mvnormal2_prior_draw_cpp(priorParams_1d, n = 3)
  expect_equal(dim(result_1d$mu), c(1, 1, 3))
  expect_equal(dim(result_1d$sig), c(1, 1, 3))

  # Test with very small nu0 (should still be valid)
  priorParams_small_nu <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2),
    phi0 = diag(2),
    nu0 = 2.1  # Just above minimum for 2D
  )

  result_small_nu <- mvnormal2_prior_draw_cpp(priorParams_small_nu, n = 2)
  expect_equal(dim(result_small_nu$mu), c(1, 2, 2))

  # Check that covariances are still positive definite
  for (i in 1:2) {
    eigenvals <- eigen(result_small_nu$sig[,,i])$values
    expect_true(all(eigenvals > 0))
  }
})

test_that("MVNormal2 cluster label changes work correctly", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Set up a simple test case
  set.seed(333)
  data <- mvtnorm::rmvnorm(10, c(0, 0), diag(2))

  dpObj <- list(
    data = data,
    n = 10,
    alpha = 1.0,
    alphaPriorParameters = c(2, 4),
    clusterLabels = c(0, 0, 0, 1, 1, 1, 2, 2, 2, 2),
    pointsPerCluster = c(3, 3, 4),
    numberClusters = 3,
    clusterParameters = list(
      mu = array(rnorm(6), dim = c(1, 2, 3)),
      sig = array(c(rep(c(1, 0, 0, 1), 3)), dim = c(2, 2, 3))
    ),
    mixingDistribution = list(
      priorParameters = list(
        mu0 = c(0, 0),
        sigma0 = diag(2),
        phi0 = diag(2),
        nu0 = 4
      )
    ),
    m = 3
  )

  # Test cluster component update WITHOUT manual modification
  # The function handles the complete update cycle internally
  result <- nonconjugate_mvnormal2_cluster_component_update_cpp(dpObj)

  # Check consistency
  expect_equal(sum(result$pointsPerCluster), 10)
  expect_true(all(result$pointsPerCluster >= 0))
  expect_equal(length(result$pointsPerCluster), result$numberClusters)

  # Check that parameters have correct dimensions
  expect_equal(dim(result$clusterParameters$mu)[3], result$numberClusters)
  expect_equal(dim(result$clusterParameters$sig)[3], result$numberClusters)
})

test_that("MVNormal2 handles numerical edge cases", {
  priorParams <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2),
    phi0 = diag(2) * 0.01,  # Very small scale
    nu0 = 4
  )

  # Generate data with very different scales
  x <- matrix(c(1e6, 1e-6, 1e6, 1e-6), ncol = 2)

  # Should not crash
  expect_error(mvnormal2_posterior_draw_cpp(priorParams, x, n = 5), NA)

  # Test with singular data (all points the same)
  x_singular <- matrix(rep(c(1, 2), 5), ncol = 2, byrow = TRUE)

  # Should handle gracefully
  result <- mvnormal2_posterior_draw_cpp(priorParams, x_singular, n = 3)
  expect_equal(dim(result$mu), c(1, 2, 3))
  expect_equal(dim(result$sig), c(2, 2, 3))
})

test_that("Hierarchical MVNormal2 integration test", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))
  skip_if_not(requireNamespace("gtools", quietly = TRUE))

  # Enable C++ samplers
  old_setting <- enable_cpp_hierarchical_samplers(TRUE)
  on.exit(enable_cpp_hierarchical_samplers(old_setting))

  # Create test data with known structure
  set.seed(444)
  true_mu1 <- c(-2, 0)
  true_mu2 <- c(2, 0)
  true_sigma <- diag(2)

  dataList <- list(
    rbind(
      mvtnorm::rmvnorm(15, true_mu1, true_sigma),
      mvtnorm::rmvnorm(10, true_mu2, true_sigma)
    ),
    rbind(
      mvtnorm::rmvnorm(20, true_mu1, true_sigma),
      mvtnorm::rmvnorm(15, true_mu2, true_sigma)
    )
  )

  # Fit hierarchical model
  hdp <- DirichletProcessHierarchicalMvnormal2(
    dataList = dataList,
    g0Priors = list(
      nu0 = 4,
      phi0 = diag(2) * 2,
      mu0 = c(0, 0),
      sigma0 = diag(2) * 3
    ),
    gammaPriors = c(2, 2),
    alphaPriors = c(2, 4),
    numSticks = 20
  )

  # Fit for a reasonable number of iterations
  hdp_fit <- Fit(hdp, its = 10, progressBar = FALSE)

  # Check basic properties
  expect_is(hdp_fit, "hierarchical")
  expect_equal(length(hdp_fit$indDP), 2)
  expect_is(hdp_fit$gammaValues, "numeric")
  expect_equal(length(hdp_fit$gammaValues), 10)

  # Check that clustering is reasonable
  # (with only 10 iterations, we don't expect perfect clustering)
  for (i in 1:2) {
    n_clusters <- hdp_fit$indDP[[i]]$numberClusters
    expect_true(n_clusters >= 1)
    expect_true(n_clusters <= nrow(dataList[[i]]))
  }
})

test_that("MVNormal2 memory safety", {
  # Test that repeated calls don't cause memory issues
  priorParams <- list(
    mu0 = c(0, 0),
    sigma0 = diag(2),
    phi0 = diag(2),
    nu0 = 4
  )

  # Repeated prior draws
  for (i in 1:10) {
    result <- mvnormal2_prior_draw_cpp(priorParams, n = 100)
    expect_equal(dim(result$mu), c(1, 2, 100))
  }

  # Repeated posterior draws
  x <- matrix(rnorm(50), ncol = 2)
  for (i in 1:10) {
    result <- mvnormal2_posterior_draw_cpp(priorParams, x, n = 50)
    expect_equal(dim(result$mu), c(1, 2, 50))
  }

  # Large hierarchical structure
  dataList <- lapply(1:10, function(i) matrix(rnorm(100), ncol = 2))

  hdp <- DirichletProcessHierarchicalMvnormal2(
    dataList = dataList,
    g0Priors = priorParams,
    numSticks = 50
  )

  # Should not crash with many DPs
  expect_equal(length(hdp$indDP), 10)
})
