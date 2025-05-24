context("Beta Distribution Full Integration Test")

test_that("Complete Beta DP workflow with C++ implementation", {
  # Enable C++ implementations
  old_setting <- enable_cpp_samplers(TRUE)

  set.seed(123)
  # Generate synthetic data with 3 Beta clusters
  maxT <- 1
  n_per_cluster <- 20

  y <- c(
    rbeta(n_per_cluster, 2, 8) * maxT,   # Mean ~0.2
    rbeta(n_per_cluster, 5, 5) * maxT,   # Mean ~0.5
    rbeta(n_per_cluster, 8, 2) * maxT    # Mean ~0.8
  )

  # Shuffle data
  y <- sample(y)

  # Create and fit DP
  dp <- DirichletProcessBeta(y, maxT, verbose = FALSE)

  # Run full MCMC
  dp <- Fit(dp, its = 50, progressBar = FALSE)

  # Check results
  expect_true(dp$numberClusters >= 2 && dp$numberClusters <= 5)
  expect_equal(sum(dp$pointsPerCluster), length(y))
  expect_true(all(dp$clusterParameters$mu > 0 & dp$clusterParameters$mu < maxT))
  expect_true(all(dp$clusterParameters$nu > 0))

  # Check convergence
  expect_true(length(unique(dp$alphaChain)) > 1)
  expect_true(all(dp$likelihoodChain > -Inf))

  # Restore setting
  enable_cpp_samplers(old_setting)
})

test_that("Performance comparison: Beta R vs C++", {
  skip_if_not(interactive(), "Performance test only run interactively")

  set.seed(456)
  y <- rbeta(100, 3, 6)

  # R implementation
  dp_r <- DirichletProcessBeta(y, 1, verbose = FALSE)
  enable_cpp_samplers(FALSE)

  time_r <- system.time({
    dp_r <- Fit(dp_r, its = 20, progressBar = FALSE)
  })

  # C++ implementation
  dp_cpp <- DirichletProcessBeta(y, 1, verbose = FALSE)
  enable_cpp_samplers(TRUE)

  time_cpp <- system.time({
    dp_cpp <- Fit(dp_cpp, its = 20, progressBar = FALSE)
  })

  cat("\nBeta Performance Comparison (20 iterations):\n")
  cat("R implementation:", round(time_r["elapsed"], 3), "seconds\n")
  cat("C++ implementation:", round(time_cpp["elapsed"], 3), "seconds\n")
  cat("Speedup:", round(time_r["elapsed"] / time_cpp["elapsed"], 1), "x\n")

  enable_cpp_samplers(FALSE)
})
