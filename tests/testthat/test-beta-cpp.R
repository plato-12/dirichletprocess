context("Beta Distribution C++ Implementation")

test_that("Beta PriorDraw C++ matches R implementation", {
  set.seed(42)
  priorParams <- c(2, 8)  # shape and rate for inverse gamma
  maxT <- 1
  n <- 10

  # Create mixing distribution object for R implementation
  mdObj <- BetaMixtureCreate(priorParams, c(1, 1), maxT)

  # R implementation
  set.seed(42)
  r_result <- PriorDraw(mdObj, n)

  # C++ implementation
  set.seed(42)
  cpp_result <- beta_prior_draw_cpp(priorParams, maxT, n)

  # Check structure
  expect_equal(names(cpp_result), c("mu", "nu"))
  expect_equal(dim(cpp_result$mu), c(1, 1, n))
  expect_equal(dim(cpp_result$nu), c(1, 1, n))

  # Check values are in valid range
  expect_true(all(cpp_result$mu > 0 & cpp_result$mu < maxT))
  expect_true(all(cpp_result$nu > 0))

  # Statistical properties should be similar
  expect_equal(mean(cpp_result$mu), mean(r_result$mu), tolerance = 0.1)
  expect_equal(mean(cpp_result$nu), mean(r_result$nu), tolerance = 0.2)
})

test_that("Beta Likelihood C++ matches R implementation", {
  priorParams <- c(2, 8)
  maxT <- 1
  mdObj <- BetaMixtureCreate(priorParams, c(1, 1), maxT)

  # Test data
  x <- seq(0.1, 0.9, by = 0.1)

  # Test parameters
  mu <- 0.5
  nu <- 4.0

  # R implementation
  theta <- list(
    mu = array(mu, dim = c(1, 1, 1)),
    nu = array(nu, dim = c(1, 1, 1))
  )
  r_lik <- Likelihood(mdObj, x, theta)

  # C++ implementation
  cpp_lik <- beta_likelihood_cpp(x, mu, nu, maxT)

  # Should match exactly
  expect_equal(cpp_lik, r_lik, tolerance = 1e-10)
})

test_that("Beta PriorDensity C++ matches R implementation", {
  priorParams <- c(2, 8)
  maxT <- 1
  mdObj <- BetaMixtureCreate(priorParams, c(1, 1), maxT)

  # Test several parameter values
  test_params <- list(
    list(mu = 0.5, nu = 2.0),
    list(mu = 0.2, nu = 10.0),
    list(mu = 0.8, nu = 0.5)
  )

  for (params in test_params) {
    # R implementation
    theta <- list(
      mu = array(params$mu, dim = c(1, 1, 1)),
      nu = array(params$nu, dim = c(1, 1, 1))
    )
    r_density <- PriorDensity(mdObj, theta)

    # C++ implementation
    cpp_density <- beta_prior_density_cpp(params$mu, params$nu, priorParams, maxT)

    # Should match closely
    expect_equal(cpp_density, r_density, tolerance = 1e-10)
  }
})

test_that("Beta PosteriorDraw C++ produces valid samples", {
  set.seed(123)
  priorParams <- c(2, 8)
  maxT <- 1
  mhStepSize <- c(0.1, 0.1)

  # Generate test data from a known Beta distribution
  a <- 3
  b <- 7
  x <- matrix(rbeta(20, a, b) * maxT, ncol = 1)

  # C++ implementation
  set.seed(456)
  cpp_result <- beta_posterior_draw_cpp(priorParams, maxT, mhStepSize, x,
                                        n = 1, mhDraws = 500)

  # Check structure
  expect_equal(names(cpp_result), c("mu", "nu"))
  expect_equal(dim(cpp_result$mu), c(1, 1, 1))
  expect_equal(dim(cpp_result$nu), c(1, 1, 1))

  # Check values are in valid range
  expect_true(all(cpp_result$mu > 0 & cpp_result$mu < maxT))
  expect_true(all(cpp_result$nu > 0))

  # The posterior mean should be close to the true mean
  true_mean <- a / (a + b) * maxT
  posterior_mu <- cpp_result$mu[1]
  expect_true(abs(posterior_mu - true_mean) < 0.2)
})

test_that("Beta Metropolis-Hastings sampler works correctly", {
  set.seed(789)
  priorParams <- c(2, 8)
  maxT <- 1
  mhStepSize <- c(0.1, 0.1)

  # Generate test data
  x <- matrix(rbeta(30, 2, 5) * maxT, ncol = 1)

  # Run MH sampler
  startMu <- 0.3
  startNu <- 5.0
  noDraws <- 1000

  mh_result <- beta_metropolis_hastings_cpp(x, startMu, startNu,
                                            priorParams, maxT,
                                            mhStepSize, noDraws)

  # Check structure
  expect_equal(names(mh_result), c("mu", "nu"))
  expect_equal(length(mh_result$mu), noDraws)
  expect_equal(length(mh_result$nu), noDraws)

  # Check convergence by looking at the second half of the chain
  mu_samples <- mh_result$mu[(noDraws/2):noDraws]
  nu_samples <- mh_result$nu[(noDraws/2):noDraws]

  # All samples should be in valid range
  expect_true(all(mu_samples > 0 & mu_samples < maxT))
  expect_true(all(nu_samples > 0))

  # Check that we have reasonable mixing (not stuck at one value)
  expect_true(length(unique(mu_samples)) > 10)
  expect_true(length(unique(nu_samples)) > 10)
})

test_that("NonconjugateBetaClusterParameterUpdate C++ works correctly", {
  set.seed(999)
  # Create a Beta DP object with known clusters
  maxT <- 1
  y1 <- rbeta(15, shape1 = 2, shape2 = 8) * maxT
  y2 <- rbeta(15, shape1 = 8, shape2 = 2) * maxT
  y <- c(y1, y2)

  dp <- DirichletProcessBeta(y, maxT, verbose = FALSE)

  # Set up known clustering
  dp$clusterLabels <- c(rep(1, 15), rep(2, 15))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(15, 15)
  dp$clusterParameters <- list(
    mu = array(c(0.3, 0.7), dim = c(1, 1, 2)),
    nu = array(c(5, 5), dim = c(1, 1, 2))
  )

  # Enable C++ samplers
  old_setting <- enable_cpp_samplers(TRUE)

  # Convert to 0-indexed for C++
  dp_cpp <- dp
  dp_cpp$clusterLabels <- dp$clusterLabels - 1

  # Run parameter update using C++
  set.seed(111)
  cpp_params <- nonconjugate_beta_cluster_parameter_update_cpp(dp_cpp)

  # Check that parameters were updated reasonably
  # Cluster 1 should have mean around 0.2 (since beta(2,8) has mean 0.2)
  expect_true(abs(cpp_params$mu[1] - 0.2) < 0.15)
  # Cluster 2 should have mean around 0.8 (since beta(8,2) has mean 0.8)
  expect_true(abs(cpp_params$mu[2] - 0.8) < 0.15)

  # Precision parameters should be positive
  expect_true(all(cpp_params$nu > 0))

  # Restore setting
  enable_cpp_samplers(old_setting)
})

test_that("End-to-end Beta C++ sampler test", {
  set.seed(2025)
  # Generate data with clear clusters
  maxT <- 1
  y <- c(rbeta(20, 2, 8) * maxT,  # Mean around 0.2
         rbeta(20, 5, 5) * maxT,  # Mean around 0.5
         rbeta(20, 8, 2) * maxT)  # Mean around 0.8

  # Initialize DP
  dp <- DirichletProcessBeta(y, maxT, verbose = FALSE)

  # Enable C++ samplers
  old_setting <- enable_cpp_samplers(TRUE)

  # Run 10 iterations
  for (iter in 1:10) {
    # Update cluster parameters using C++
    dp <- ClusterParameterUpdate(dp)

    # Basic sanity checks
    expect_true(all(dp$clusterParameters$mu > 0 & dp$clusterParameters$mu < maxT))
    expect_true(all(dp$clusterParameters$nu > 0))
  }

  # Check that we have reasonable parameter values
  final_mus <- as.numeric(dp$clusterParameters$mu)

  cat("\nEnd-to-end Beta test summary:\n")
  cat("Number of clusters:", dp$numberClusters, "\n")
  cat("Cluster sizes:", as.numeric(dp$pointsPerCluster), "\n")
  cat("Cluster means (mu):", round(final_mus, 3), "\n")
  cat("Cluster precisions (nu):", round(as.numeric(dp$clusterParameters$nu), 3), "\n")

  # Restore setting
  enable_cpp_samplers(old_setting)
})

test_that("C++ and R implementations produce similar results", {
  set.seed(333)
  maxT <- 1
  priorParams <- c(2, 8)
  mhStepSize <- c(0.1, 0.1)

  # Generate test data
  x <- matrix(rbeta(25, 3, 6) * maxT, ncol = 1)

  # Create mixing distribution
  mdObj <- BetaMixtureCreate(priorParams, mhStepSize, maxT)

  # R implementation
  set.seed(444)
  r_result <- PosteriorDraw(mdObj, x, n = 1)

  # C++ implementation
  set.seed(444)
  cpp_result <- beta_posterior_draw_cpp(priorParams, maxT, mhStepSize, x,
                                        n = 1, mhDraws = 250)

  # Results won't be identical due to RNG differences, but should be similar
  # Check that both are drawing from similar regions of parameter space
  r_mu <- r_result$mu[1]
  cpp_mu <- cpp_result$mu[1]

  r_nu <- r_result$nu[1]
  cpp_nu <- cpp_result$nu[1]

  # They should be in the same ballpark
  expect_true(abs(r_mu - cpp_mu) < 0.3)
  expect_true(abs(log(r_nu) - log(cpp_nu)) < 1.0)  # Use log scale for nu
})
