context("Beta DP Performance")

test_that("C++ implementation is faster than R", {
  skip_if_not(can_use_cpp(DirichletProcessBeta(rbeta(10, 2, 2))),
              "C++ Beta not available")
  skip_if_not(requireNamespace("microbenchmark", quietly = TRUE))

  set.seed(4567)
  data <- rbeta(100, 3, 7)

  # Time R implementation
  set_use_cpp(FALSE)
  time_r <- system.time({
    dp_r <- DirichletProcessBeta(data, verbose = FALSE)
    dp_r <- Fit(dp_r, its = 50, progressBar = FALSE)
  })["elapsed"]

  # Time C++ implementation
  set_use_cpp(TRUE)
  time_cpp <- system.time({
    dp_cpp <- DirichletProcessBeta(data, verbose = FALSE)
    dp_cpp <- Fit(dp_cpp, its = 50, progressBar = FALSE)
  })["elapsed"]

  # C++ should be at least 5x faster
  speedup <- time_r / time_cpp
  expect_true(speedup > 5,
              info = sprintf("Speedup: %.1fx (R: %.3fs, C++: %.3fs)",
                             speedup, time_r, time_cpp))
})

test_that("Beta DP scales well with data size", {
  skip_on_cran()  # Skip on CRAN due to time

  set_use_cpp(TRUE)

  # Test scaling
  sizes <- c(50, 100, 200, 400)
  times <- numeric(length(sizes))

  for (i in seq_along(sizes)) {
    n <- sizes[i]
    data <- rbeta(n, 2, 8)

    times[i] <- system.time({
      dp <- DirichletProcessBeta(data, verbose = FALSE)
      dp <- Fit(dp, its = 20, progressBar = FALSE)
    })["elapsed"]
  }

  # Check approximately linear scaling
  # Time should roughly double when data doubles
  scaling_factors <- times[-1] / times[-length(times)]
  expect_true(all(scaling_factors < 3))  # Less than cubic scaling
})
