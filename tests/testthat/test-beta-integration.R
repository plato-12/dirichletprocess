context("Beta Distribution Full Integration Test")

test_that("Complete Beta DP workflow (expecting R fallback for C++)", {
  old_setting <- enable_cpp_samplers(TRUE) # Enable C++ path (will hit stubs)
  set.seed(123)
  maxT <- 1
  n_per_cluster <- 20
  y <- c(
    rbeta(n_per_cluster, 2, 8) * maxT,
    rbeta(n_per_cluster, 5, 5) * maxT,
    rbeta(n_per_cluster, 8, 2) * maxT
  )
  y <- sample(y)

  # Create DP with mhDraws specified for non-conjugate samplers
  dp <- DirichletProcessBeta(y, maxT, verbose = FALSE, mhDraws = 50)

  # Expect warnings as C++ stubs for non-conjugate parts are called
  # and R fallbacks are used.
  expect_warning(
    dp_fitted <- Fit(dp, its = 20, progressBar = FALSE), # Reduced iterations for faster test
    regexp = "C\\+\\+ function 'nonconjugate_beta_cluster_parameter_update_cpp' is a STUB"
  )
  # If ClusterComponentUpdate also has a C++ stub for beta:
  # expect_warning(
  #   dp_fitted <- Fit(dp, its = 20, progressBar = FALSE),
  #   regexp = "C\\+\\+ function 'nonconjugate_beta_cluster_component_update_cpp' is a STUB" # Or similar
  # )


  expect_true(dp_fitted$numberClusters >= 1 && dp_fitted$numberClusters <= (length(y)/2)) # More relaxed check
  expect_equal(sum(dp_fitted$pointsPerCluster), length(y))
  if (dp_fitted$numberClusters > 0) {
    expect_true(all(dp_fitted$clusterParameters$mu > 0 & dp_fitted$clusterParameters$mu < maxT))
    expect_true(all(dp_fitted$clusterParameters$nu > 0))
  }
  expect_true(length(unique(dp_fitted$alphaChain)) > 1 || length(dp_fitted$alphaChain) == 1) # Can be 1 if alpha fixed
  expect_true(all(is.finite(dp_fitted$likelihoodChain))) # Check for finite likelihood

  enable_cpp_samplers(old_setting) # Restore original setting
})

test_that("Performance comparison: Beta R vs C++ (Stub)", {
  skip_if_not(interactive(), "Performance test only run interactively")
  set.seed(456)
  y <- rbeta(50, 3, 6) # Reduced data size for faster interactive test

  # R implementation
  dp_r <- DirichletProcessBeta(y, 1, verbose = FALSE, mhDraws = 30)
  enable_cpp_samplers(FALSE)
  time_r <- system.time({
    dp_r <- Fit(dp_r, its = 10, progressBar = FALSE) # Reduced iterations
  })

  # C++ implementation (will use R fallback due to stubs)
  dp_cpp <- DirichletProcessBeta(y, 1, verbose = FALSE, mhDraws = 30)
  enable_cpp_samplers(TRUE)
  # Expect warnings from the C++ stubs
  time_cpp <- system.time({
    expect_warning(dp_cpp <- Fit(dp_cpp, its = 10, progressBar = FALSE))
  })

  cat("\nBeta Performance Comparison (10 iterations, 50 data points - C++ uses R fallback):\n")
  cat("R implementation:", round(time_r["elapsed"], 3), "seconds\n")
  cat("C++ (fallback) implementation:", round(time_cpp["elapsed"], 3), "seconds\n")
  # Speedup might not be meaningful here as C++ uses R fallback.
  cat("Speedup (R vs C++ fallback):", round(time_r["elapsed"] / time_cpp["elapsed"], 1), "x\n")

  enable_cpp_samplers(FALSE) # Important to reset
})
