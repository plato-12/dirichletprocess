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
