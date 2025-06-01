test_that("Hierarchical Beta DP C++ implementation works", {
  skip_if_not_installed("gtools")
  skip_if_not_installed("mvtnorm")

  # Save current state and ensure clean start
  old_state <- using_cpp_hierarchical_samplers()
  on.exit(enable_cpp_hierarchical_samplers(old_state), add = TRUE)
  enable_cpp_hierarchical_samplers(FALSE)

  # Generate small test data
  set.seed(123)
  dataList <- list(
    rbeta(10, 2, 5),
    rbeta(10, 5, 2)
  )

  # Test that C++ implementation can be enabled/disabled
  expect_false(using_cpp_hierarchical_samplers())
  enable_cpp_hierarchical_samplers(TRUE)
  expect_true(using_cpp_hierarchical_samplers())
  enable_cpp_hierarchical_samplers(FALSE)
  expect_false(using_cpp_hierarchical_samplers())
})

test_that("Hierarchical Beta DP C++ and R implementations produce valid results", {
  skip_if_not_installed("gtools")
  skip_if_not_installed("mvtnorm")

  # Ensure clean state
  old_state <- using_cpp_hierarchical_samplers()
  on.exit(enable_cpp_hierarchical_samplers(old_state), add = TRUE)
  enable_cpp_hierarchical_samplers(FALSE)

  set.seed(42)
  dataList <- list(
    rbeta(15, 2, 5),
    rbeta(15, 5, 2)
  )

  # Create DP with R implementation
  dp_r <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 10,
    mhDraws = 10
  )

  # Fit with R implementation
  set.seed(100)
  dp_r_fit <- Fit(dp_r, its = 2, progressBar = FALSE)  # Reduced iterations for testing

  # Enable C++ and create DP
  enable_cpp_hierarchical_samplers(TRUE)

  # Test that creation works
  expect_error({
    dp_cpp <- DirichletProcessHierarchicalBeta(
      dataList = dataList,
      maxY = 1,
      priorParameters = c(2, 8),
      hyperPriorParameters = c(1, 0.125),
      gammaPriors = c(2, 4),
      alphaPriors = c(2, 4),
      mhStepSize = c(0.1, 0.1),
      numSticks = 10,
      mhDraws = 10
    )
  }, NA)  # Expect no error

  # Only proceed with fitting if C++ implementation is available
  # Check if the C++ functions exist before testing
  if (exists("hierarchical_beta_fit_cpp")) {
    dp_cpp <- DirichletProcessHierarchicalBeta(
      dataList = dataList,
      maxY = 1,
      priorParameters = c(2, 8),
      hyperPriorParameters = c(1, 0.125),
      gammaPriors = c(2, 4),
      alphaPriors = c(2, 4),
      mhStepSize = c(0.1, 0.1),
      numSticks = 10,
      mhDraws = 10
    )

    # Fit with C++ implementation
    set.seed(100)
    expect_error({
      dp_cpp_fit <- Fit(dp_cpp, its = 2, progressBar = FALSE)
    }, NA)  # Expect no error

    # Basic structure tests
    expect_equal(length(dp_r_fit$indDP), length(dp_cpp_fit$indDP))
    expect_equal(length(dp_r_fit$globalParameters), length(dp_cpp_fit$globalParameters))

    # Check that gamma values are numeric
    expect_true(is.numeric(dp_r_fit$gamma))
    expect_true(is.numeric(dp_cpp_fit$gamma))

    # Check that both have valid cluster assignments
    for (i in seq_along(dataList)) {
      expect_equal(length(dp_r_fit$indDP[[i]]$clusterLabels),
                   length(dataList[[i]]))
      expect_equal(length(dp_cpp_fit$indDP[[i]]$clusterLabels),
                   length(dataList[[i]]))
      expect_true(all(dp_r_fit$indDP[[i]]$clusterLabels >= 1))
      expect_true(all(dp_cpp_fit$indDP[[i]]$clusterLabels >= 1))
    }
  } else {
    skip("C++ implementation not available")
  }
})

test_that("Individual update functions work with C++ implementation", {
  skip_if_not_installed("gtools")
  skip_if_not_installed("mvtnorm")

  # Ensure clean state
  old_state <- using_cpp_hierarchical_samplers()
  on.exit(enable_cpp_hierarchical_samplers(old_state), add = TRUE)

  # Skip if C++ functions not available
  if (!exists("hierarchical_beta_cluster_component_update_cpp")) {
    skip("C++ implementation not available")
  }

  set.seed(123)
  dataList <- list(
    rbeta(10, 2, 5),
    rbeta(10, 5, 2)
  )

  enable_cpp_hierarchical_samplers(TRUE)

  dp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 10,
    mhDraws = 5
  )

  # Test ClusterComponentUpdate
  expect_error({
    dp_updated <- ClusterComponentUpdate(dp)
  }, NA)
  expect_s3_class(dp_updated, "hierarchical")
  expect_equal(length(dp_updated$indDP), length(dataList))

  # Test GlobalParameterUpdate
  expect_error({
    dp_updated <- GlobalParameterUpdate(dp)
  }, NA)
  expect_s3_class(dp_updated, "hierarchical")
  expect_true(!is.null(dp_updated$globalParameters))

  # Test UpdateG0
  expect_error({
    dp_updated <- UpdateG0(dp)
  }, NA)
  expect_s3_class(dp_updated, "hierarchical")
  expect_true(!is.null(dp_updated$globalStick))

  # Test UpdateGamma
  expect_error({
    dp_updated <- UpdateGamma(dp)
  }, NA)
  expect_s3_class(dp_updated, "hierarchical")
  expect_true(is.numeric(dp_updated$gamma))
  expect_true(dp_updated$gamma > 0)
})

test_that("C++ implementation handles edge cases", {
  skip_if_not_installed("gtools")
  skip_if_not_installed("mvtnorm")

  # Ensure clean state
  old_state <- using_cpp_hierarchical_samplers()
  on.exit(enable_cpp_hierarchical_samplers(old_state), add = TRUE)

  # Skip if C++ functions not available
  if (!exists("hierarchical_beta_fit_cpp")) {
    skip("C++ implementation not available")
  }

  enable_cpp_hierarchical_samplers(TRUE)

  # Test with very small dataset
  dataList <- list(
    rbeta(3, 2, 5),
    rbeta(3, 5, 2)
  )

  dp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 5,
    mhDraws = 5
  )

  # Should not crash
  expect_error(Fit(dp, its = 2, progressBar = FALSE), NA)

  # Test with single dataset
  dataList_single <- list(rbeta(10, 2, 5))

  dp_single <- DirichletProcessHierarchicalBeta(
    dataList = dataList_single,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 5,
    mhDraws = 5
  )

  expect_error(Fit(dp_single, its = 2, progressBar = FALSE), NA)
})

test_that("C++ mixing distribution creation works", {
  skip_if_not_installed("gtools")
  skip_if_not_installed("mvtnorm")

  # Ensure clean state
  old_state <- using_cpp_hierarchical_samplers()
  on.exit(enable_cpp_hierarchical_samplers(old_state), add = TRUE)

  # Skip if C++ functions not available
  if (!exists("hierarchical_beta_mixing_create_cpp")) {
    skip("C++ implementation not available")
  }

  enable_cpp_hierarchical_samplers(TRUE)

  # Test HierarchicalBetaCreate
  expect_error({
    mdobj_list <- HierarchicalBetaCreate(
      n = 2,
      priorParameters = c(2, 8),
      hyperPriorParameters = c(1, 0.125),
      alphaPrior = c(2, 4),
      maxT = 1,
      gammaPrior = c(2, 4),
      mhStepSize = c(0.1, 0.1),
      num_sticks = 10
    )
  }, NA)

  expect_equal(length(mdobj_list), 2)
  expect_true(all(sapply(mdobj_list, function(x) "hierarchical" %in% class(x))))
  expect_true(all(sapply(mdobj_list, function(x) "beta" %in% class(x))))
  expect_true(all(sapply(mdobj_list, function(x) !is.null(x$theta_k))))
  expect_true(all(sapply(mdobj_list, function(x) !is.null(x$beta_k))))
  expect_true(all(sapply(mdobj_list, function(x) !is.null(x$pi_k))))
})

test_that("Hierarchical Beta DP recovers known parameters", {
  skip_if_not_installed("gtools")

  # Generate data from known clusters
  set.seed(123)
  true_mu1 <- 0.3
  true_mu2 <- 0.7
  true_nu <- 10

  # Generate data from two clear clusters
  group1_cluster1 <- rbeta(50, true_mu1 * true_nu, (1 - true_mu1) * true_nu)
  group1_cluster2 <- rbeta(50, true_mu2 * true_nu, (1 - true_mu2) * true_nu)

  group2_cluster1 <- rbeta(40, true_mu1 * true_nu, (1 - true_mu1) * true_nu)
  group2_cluster2 <- rbeta(60, true_mu2 * true_nu, (1 - true_mu2) * true_nu)

  dataList <- list(
    c(group1_cluster1, group1_cluster2),
    c(group2_cluster1, group2_cluster2)
  )

  enable_cpp_hierarchical_samplers(TRUE)

  dp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 50
  )

  # Run for more iterations
  dp_fit <- Fit(dp, its = 1000, progressBar = FALSE)

  # Check that we found approximately 2 global clusters
  global_params <- dp_fit$globalParameters
  mu_values <- global_params[[1]]

  # Find unique clusters (accounting for MCMC variability)
  unique_mus <- unique(round(mu_values, 1))
  expect_true(length(unique_mus) >= 2 && length(unique_mus) <= 4)

  # Check that recovered parameters are close to true values
  expect_true(any(abs(mu_values - true_mu1) < 0.1))
  expect_true(any(abs(mu_values - true_mu2) < 0.1))
})

test_that("MCMC chains show convergence", {
  skip_if_not_installed("gtools")

  set.seed(42)
  dataList <- list(
    rbeta(50, 2, 5),
    rbeta(50, 5, 2)
  )

  enable_cpp_hierarchical_samplers(TRUE)

  dp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 20
  )

  dp_fit <- Fit(dp, its = 500, progressBar = FALSE)

  # Check gamma chain convergence
  gamma_chain <- dp_fit$gammaValues

  # Simple convergence check: second half should have similar mean to first half
  first_half_mean <- mean(gamma_chain[1:250])
  second_half_mean <- mean(gamma_chain[251:500])

  expect_true(abs(first_half_mean - second_half_mean) / first_half_mean < 0.2)

  # Check that gamma values are reasonable
  expect_true(all(gamma_chain > 0))
  expect_true(all(gamma_chain < 100))  # Reasonable upper bound
})

test_that("Prior parameters update correctly", {
  skip_if_not_installed("gtools")

  set.seed(123)
  dataList <- list(
    rbeta(30, 2, 8),  # Data matching the prior
    rbeta(30, 2, 8)
  )

  enable_cpp_hierarchical_samplers(TRUE)

  dp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(1, 1),  # Weak prior
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 20
  )

  # Fit with prior updates
  dp_fit <- Fit(dp, its = 100, updatePrior = TRUE, progressBar = FALSE)

  # Check that prior parameters have been updated
  final_prior <- dp_fit$indDP[[1]]$mixingDistribution$priorParameters
  initial_prior <- c(1, 1)

  expect_false(all(final_prior == initial_prior))

  # The beta parameter for nu should have increased from the data
  expect_true(final_prior[2] > initial_prior[2])
})

test_that("R and C++ implementations produce statistically similar results", {
  skip_if_not_installed("gtools")

  set.seed(999)
  dataList <- list(
    rbeta(100, 3, 7),
    rbeta(100, 7, 3)
  )

  # Run with R implementation
  enable_cpp_hierarchical_samplers(FALSE)
  dp_r <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 30,
    mhDraws = 100
  )

  set.seed(100)
  dp_r_fit <- Fit(dp_r, its = 200, progressBar = FALSE)

  # Run with C++ implementation
  enable_cpp_hierarchical_samplers(TRUE)
  dp_cpp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 30,
    mhDraws = 100
  )

  set.seed(100)
  dp_cpp_fit <- Fit(dp_cpp, its = 200, progressBar = FALSE)

  # Compare posterior means of gamma
  r_gamma_mean <- mean(dp_r_fit$gammaValues[101:200])
  cpp_gamma_mean <- mean(dp_cpp_fit$gammaValues[101:200])

  # Should be within 20% of each other
  expect_true(abs(r_gamma_mean - cpp_gamma_mean) / r_gamma_mean < 0.2)

  # Compare number of clusters
  r_clusters <- sapply(dp_r_fit$indDP, function(x) x$numberClusters)
  cpp_clusters <- sapply(dp_cpp_fit$indDP, function(x) x$numberClusters)

  expect_equal(mean(r_clusters), mean(cpp_clusters), tolerance = 1)
})

test_that("C++ implementation is faster than R", {
  skip_if_not_installed("gtools")
  skip_if_not_installed("microbenchmark")

  dataList <- list(
    rbeta(200, 2, 5),
    rbeta(200, 5, 2),
    rbeta(200, 3, 3)
  )

  # Create DP object once
  dp <- DirichletProcessHierarchicalBeta(
    dataList = dataList,
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 20,
    mhDraws = 50
  )

  # Time R implementation
  enable_cpp_hierarchical_samplers(FALSE)
  r_time <- system.time({
    Fit(dp, its = 50, progressBar = FALSE)
  })["elapsed"]

  # Time C++ implementation
  enable_cpp_hierarchical_samplers(TRUE)
  cpp_time <- system.time({
    Fit(dp, its = 50, progressBar = FALSE)
  })["elapsed"]

  # C++ should be faster
  expect_true(cpp_time < r_time)

  # Print speedup for information
  cat("\nSpeedup:", round(r_time / cpp_time, 2), "x\n")
})

test_that("Implementation handles numerical edge cases", {
  skip_if_not_installed("gtools")

  enable_cpp_hierarchical_samplers(TRUE)

  # Test with data very close to boundaries
  dataList <- list(
    c(rep(0.001, 10), rep(0.999, 10)),
    c(rep(0.001, 10), rep(0.999, 10))
  )

  expect_error({
    dp <- DirichletProcessHierarchicalBeta(
      dataList = dataList,
      maxY = 1,
      priorParameters = c(2, 8),
      hyperPriorParameters = c(1, 0.125),
      gammaPriors = c(2, 4),
      alphaPriors = c(2, 4),
      mhStepSize = c(0.01, 0.01),  # Small step size for stability
      numSticks = 10
    )

    dp_fit <- Fit(dp, its = 20, progressBar = FALSE)
  }, NA)  # Should not error

  # Check outputs are finite
  expect_true(all(is.finite(dp_fit$gammaValues)))

  for (i in seq_along(dp_fit$indDP)) {
    params <- dp_fit$indDP[[i]]$clusterParameters
    expect_true(all(is.finite(params[[1]])))
    expect_true(all(is.finite(params[[2]])))
  }
})

test_that("No memory leaks in repeated fitting", {
  skip_if_not_installed("gtools")
  skip_on_cran()  # Memory tests can be flaky on CRAN

  enable_cpp_hierarchical_samplers(TRUE)

  dataList <- list(
    rbeta(100, 2, 5),
    rbeta(100, 5, 2)
  )

  # Run garbage collection first
  gc()

  # Get initial memory
  mem_before <- as.numeric(gc()[2, 2])

  # Run multiple fits
  for (i in 1:10) {
    dp <- DirichletProcessHierarchicalBeta(
      dataList = dataList,
      maxY = 1,
      priorParameters = c(2, 8),
      hyperPriorParameters = c(1, 0.125),
      gammaPriors = c(2, 4),
      alphaPriors = c(2, 4),
      mhStepSize = c(0.1, 0.1),
      numSticks = 20
    )

    dp_fit <- Fit(dp, its = 50, progressBar = FALSE)
    rm(dp, dp_fit)
  }

  # Force garbage collection
  gc()

  # Check memory after
  mem_after <- as.numeric(gc()[2, 2])

  # Memory increase should be minimal (less than 50MB)
  expect_true((mem_after - mem_before) < 50)
})








































































