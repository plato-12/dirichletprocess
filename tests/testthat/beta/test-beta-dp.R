context("Beta Dirichlet Process")

test_that("Beta DP initialization works correctly", {
  set.seed(100)
  data <- rbeta(20, 2, 8)

  dp <- DirichletProcessBeta(
    y = data,
    maxY = 1,
    g0Priors = c(2, 8),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    verbose = FALSE
  )

  expect_s3_class(dp, "dirichletprocess")
  expect_s3_class(dp, "beta")
  expect_s3_class(dp, "nonconjugate")

  expect_equal(dp$n, 20)
  expect_equal(dp$data, matrix(data, ncol = 1))
  expect_equal(dp$numberClusters, 1)
  expect_equal(dp$clusterLabels, rep(1, 20))
  expect_equal(dp$pointsPerCluster, 20)

  # Check parameters are properly initialized
  expect_length(dp$clusterParameters, 2)
  expect_true(all(dp$clusterParameters[[1]] > 0))
  expect_true(all(dp$clusterParameters[[1]] < 1))
  expect_true(all(dp$clusterParameters[[2]] > 0))
})

test_that("Beta DP cluster updates work correctly", {
  set.seed(200)
  # Create data with two clear clusters
  data <- c(rbeta(30, 2, 8), rbeta(30, 8, 2))

  dp <- DirichletProcessBeta(data, verbose = FALSE)

  # Test cluster label update
  dp_updated <- ClusterComponentUpdate(dp)

  expect_equal(length(dp_updated$clusterLabels), 60)
  expect_equal(sum(dp_updated$pointsPerCluster), 60)
  expect_true(dp_updated$numberClusters >= 1)

  # Run multiple updates to allow clustering
  for (i in 1:10) {
    dp_updated <- ClusterComponentUpdate(dp_updated)
  }

  # Should find evidence of multiple clusters
  expect_true(dp_updated$numberClusters > 1)
})

test_that("Beta DP parameter updates preserve validity", {
  set.seed(300)
  data <- rbeta(40, 3, 3)

  dp <- DirichletProcessBeta(data, verbose = FALSE)

  # Initialize with 2 clusters
  dp$clusterLabels <- c(rep(1, 20), rep(2, 20))
  dp$numberClusters <- 2
  dp$pointsPerCluster <- c(20, 20)
  dp$clusterParameters <- PriorDraw(dp$mixingDistribution, 2)

  # Update parameters
  dp_updated <- ClusterParameterUpdate(dp)

  # Check validity
  expect_length(dp_updated$clusterParameters[[1]], 2)
  expect_length(dp_updated$clusterParameters[[2]], 2)

  # All mu values should be in (0, maxT)
  mu_vals <- dp_updated$clusterParameters[[1]]
  expect_true(all(mu_vals > 0 & mu_vals < dp$mixingDistribution$maxT))

  # All nu values should be positive
  nu_vals <- dp_updated$clusterParameters[[2]]
  expect_true(all(nu_vals > 0))
})
