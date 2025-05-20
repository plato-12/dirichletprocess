context("C++ function testing")

test_that("Normal likelihood C++ implementation works correctly", {
  # Create test data
  x <- seq(-3, 3, by = 0.1)
  mdObj <- list()  # Simplified mixing object for testing
  theta <- list(0, 1)  # Mean 0, standard deviation 1

  # Get results from C++ implementation
  cpp_result <- likelihood_normal_cpp(mdObj, x, theta)

  # Calculate expected values directly in R
  r_result <- dnorm(x, 0, 1)

  # Verify results match
  expect_equal(cpp_result, r_result, tolerance = 1e-12)
})
