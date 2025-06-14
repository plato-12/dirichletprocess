# tests/testthat/test-cpp-mcmc.R

test_that("C++ and R implementations give similar results for Gaussian", {
  set.seed(123)

  # Generate test data
  true_clusters <- c(rep(1, 50), rep(2, 50))
  data <- c(rnorm(50, -2, 1), rnorm(50, 2, 1))

  # Fit with R implementation
  dp_r <- DirichletProcessGaussian(data)
  set_use_cpp(FALSE)
  dp_r <- Fit(dp_r, 1000, 500, 1)

  # Fit with C++ implementation
  dp_cpp <- DirichletProcessGaussian(data)
  set_use_cpp(TRUE)
  dp_cpp <- Fit(dp_cpp, 1000, 500, 1)

  # Check that results are similar (not identical due to randomness)
  expect_equal(length(unique(dp_r$cluster_labels[[500]])),
               length(unique(dp_cpp$cluster_labels[[500]])),
               tolerance = 1)

  expect_equal(mean(dp_r$alpha), mean(dp_cpp$alpha), tolerance = 0.5)
})

test_that("C++ implementation handles edge cases", {
  # Single observation
  dp <- DirichletProcessGaussian(c(1.0))
  set_use_cpp(TRUE)
  expect_no_error(Fit(dp, 100))

  # Empty clusters are properly removed
  # Extreme alpha values
  # etc.
})
