test_that("Hierarchical Beta DP C++ implementation works", {
  skip_if_not_installed("gtools")

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

test_that("Hierarchical Beta DP C++ and R implementations produce similar results", {
  skip_if_not_installed("gtools")

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
  dp_r_fit <- Fit(dp_r, its = 5, progressBar = FALSE)

  # Enable C++ and create DP
  enable_cpp_hierarchical_samplers(TRUE)

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
  dp_cpp_fit <- Fit(dp_cpp, its = 5, progressBar = FALSE)

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

  # Disable C++ implementations
  enable_cpp_hierarchical_samplers(FALSE)
})

test_that("Individual update functions work with C++ implementation", {
  skip_if_not_installed("gtools")

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
  dp_updated <- ClusterComponentUpdate(dp)
  expect_s3_class(dp_updated, "hierarchical")
  expect_equal(length(dp_updated$indDP), length(dataList))

  # Test GlobalParameterUpdate
  dp_updated <- GlobalParameterUpdate(dp)
  expect_s3_class(dp_updated, "hierarchical")
  expect_true(!is.null(dp_updated$globalParameters))

  # Test UpdateG0
  dp_updated <- UpdateG0(dp)
  expect_s3_class(dp_updated, "hierarchical")
  expect_true(!is.null(dp_updated$globalStick))

  # Test UpdateGamma
  dp_updated <- UpdateGamma(dp)
  expect_s3_class(dp_updated, "hierarchical")
  expect_true(is.numeric(dp_updated$gamma))
  expect_true(dp_updated$gamma > 0)

  enable_cpp_hierarchical_samplers(FALSE)
})

test_that("C++ implementation handles edge cases", {
  skip_if_not_installed("gtools")

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

  enable_cpp_hierarchical_samplers(FALSE)
})

test_that("C++ mixing distribution creation works", {
  skip_if_not_installed("gtools")

  enable_cpp_hierarchical_samplers(TRUE)

  # Test HierarchicalBetaCreate
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

  expect_equal(length(mdobj_list), 2)
  expect_true(all(sapply(mdobj_list, function(x) "hierarchical" %in% class(x))))
  expect_true(all(sapply(mdobj_list, function(x) "beta" %in% class(x))))
  expect_true(all(sapply(mdobj_list, function(x) !is.null(x$theta_k))))
  expect_true(all(sapply(mdobj_list, function(x) !is.null(x$beta_k))))
  expect_true(all(sapply(mdobj_list, function(x) !is.null(x$pi_k))))

  enable_cpp_hierarchical_samplers(FALSE)
})
