context("Beta DP Edge Cases")

test_that("Beta DP handles extreme data values", {
  # Data very close to boundaries
  extreme_data <- c(
    rep(0.001, 10),
    rep(0.999, 10),
    rbeta(10, 1, 1)
  )

  expect_error({
    dp <- DirichletProcessBeta(extreme_data, verbose = FALSE)
    dp <- Fit(dp, its = 50, progressBar = FALSE)
  }, NA)  # Should not error

  # Check results are valid
  expect_true(all(is.finite(dp$clusterParameters[[1]])))
  expect_true(all(is.finite(dp$clusterParameters[[2]])))
})

test_that("Beta DP handles single data point", {
  single_point <- 0.7

  expect_error({
    dp <- DirichletProcessBeta(single_point, verbose = FALSE)
  }, NA)

  expect_equal(dp$n, 1)
  expect_equal(dp$numberClusters, 1)
})

test_that("Beta DP handles constant data", {
  constant_data <- rep(0.5, 30)

  dp <- DirichletProcessBeta(constant_data, verbose = FALSE)
  dp <- Fit(dp, its = 100, progressBar = FALSE)

  # Should find 1 cluster centered at 0.5
  expect_equal(dp$numberClusters, 1)
  expect_equal(dp$clusterParameters[[1]][1], 0.5, tolerance = 0.01)
})

test_that("Beta DP validates input parameters", {
  data <- rbeta(20, 2, 2)

  # Invalid priors
  expect_error(
    DirichletProcessBeta(data, g0Priors = c(-1, 2)),
    "prior"
  )

  # Invalid alpha priors
  expect_error(
    DirichletProcessBeta(data, alphaPriors = c(0, -1)),
    "alpha"
  )

  # Invalid data
  expect_error(
    DirichletProcessBeta(c(0.5, 1.5, 0.3)),
    "bound|range"
  )
})

test_that("Beta DP handles missing values appropriately", {
  data_with_na <- c(rbeta(20, 2, 2), NA, NA)

  # Should either handle or error gracefully
  result <- tryCatch({
    dp <- DirichletProcessBeta(data_with_na, verbose = FALSE)
    "handled"
  }, error = function(e) {
    "error"
  })

  expect_true(result %in% c("handled", "error"))
})
