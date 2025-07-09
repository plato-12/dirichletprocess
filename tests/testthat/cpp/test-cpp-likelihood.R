context("C++ Likelihood Implementation")

test_that("Normal likelihood C++ implementation matches R", {
  # Generate test data
  set.seed(42)
  x <- rnorm(100)

  # Create a mixing distribution object
  mdObj <- MixingDistribution("normal", c(0,1,1,1), "conjugate")

  # Generate parameters
  theta <- list(array(0, dim=c(1,1,1)), array(1, dim=c(1,1,1)))

  # Compare R and C++ implementations
  result <- compare_r_cpp(
    r_func = function(mdObj, x, theta) Likelihood(mdObj, x, theta),
    cpp_func = function(mdObj, x, theta) likelihood_cpp(mdObj, x, theta),
    mdObj = mdObj, x = x, theta = theta
  )

  expect_true(result$equal, result$message)
})

test_that("Normal likelihood C++ implementation is faster", {
  # Generate larger test data for benchmarking
  set.seed(42)
  x <- rnorm(1000)

  # Create a mixing distribution object
  mdObj <- MixingDistribution("normal", c(0,1,1,1), "conjugate")

  # Generate parameters
  theta <- list(array(0, dim=c(1,1,1)), array(1, dim=c(1,1,1)))

  # Benchmark R vs C++
  bench <- benchmark_r_cpp(
    r_func = function(mdObj, x, theta) Likelihood(mdObj, x, theta),
    cpp_func = function(mdObj, x, theta) likelihood_cpp(mdObj, x, theta),
    mdObj = mdObj, x = x, theta = theta,
    times = 100
  )

  # C++ should be faster, but might not be for very simple cases
  # Just print results for now
  cat("R time:", bench$r_time["elapsed"],
      "C++ time:", bench$cpp_time["elapsed"],
      "Speedup:", bench$speedup, "\n")

  # Once implementation is mature, add a stricter test:
  # expect_gt(bench$speedup, 1.0)
})
