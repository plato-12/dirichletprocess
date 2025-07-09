# inst/benchmarks/test-exponential-benchmark.R
# Test script for Exponential Distribution Performance
# This demonstrates the comprehensive benchmarking and memory profiling

library(dirichletprocess)
library(testthat)

# Test the memory profiling functionality
test_that("C++ memory profiling works correctly", {
  skip_if_not(using_cpp(), "C++ implementation not enabled")

  # Clear any previous tracking
  clear_memory_tracking()

  # Create some data
  set.seed(123)
  y <- rexp(100, rate = 2)

  # Fit with C++ implementation
  set_use_cpp(TRUE)
  dp <- DirichletProcessExponential(y)
  dp <- Fit(dp, 50, progressBar = FALSE)

  # Get memory tracking
  mem_track <- get_memory_tracking()

  # Check that we have memory tracking data
  expect_s3_class(mem_track, "data.frame")
  if (nrow(mem_track) > 0) {
    expect_true("description" %in% names(mem_track))
    expect_true("bytes" %in% names(mem_track))
    expect_true("mb" %in% names(mem_track))
  }
})

# Quick benchmark function
quick_exponential_benchmark <- function() {
  cat("\n=== Quick Exponential Distribution Benchmark ===\n\n")

  # Test parameters
  n_obs <- c(100, 500, 1000)
  n_iter <- 100
  n_reps <- 3

  results <- data.frame()

  for (n in n_obs) {
    cat(sprintf("Testing n = %d:\n", n))

    # Generate test data
    set.seed(42)
    data <- c(rexp(n/2, rate = 2), rexp(n/2, rate = 5))

    # Benchmark R implementation
    times_r <- numeric(n_reps)
    for (i in 1:n_reps) {
      set_use_cpp(FALSE)
      time_r <- system.time({
        dp_r <- DirichletProcessExponential(data)
        dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
      })
      times_r[i] <- time_r["elapsed"]
    }

    # Benchmark C++ implementation
    times_cpp <- numeric(n_reps)
    mem_usage <- list()

    for (i in 1:n_reps) {
      set_use_cpp(TRUE)
      clear_memory_tracking()

      time_cpp <- system.time({
        dp_cpp <- DirichletProcessExponential(data)
        dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
      })
      times_cpp[i] <- time_cpp["elapsed"]
      mem_usage[[i]] <- get_memory_tracking()
    }

    # Calculate statistics
    mean_r <- mean(times_r)
    mean_cpp <- mean(times_cpp)
    speedup <- mean_r / mean_cpp

    # Memory usage from last run
    total_mem_mb <- 0
    if (length(mem_usage[[n_reps]]) > 0 && nrow(mem_usage[[n_reps]]) > 0) {
      total_mem_mb <- sum(mem_usage[[n_reps]]$mb)
    }

    cat(sprintf("  R: %.3fs, C++: %.3fs, Speedup: %.2fx, Memory: %.2f MB\n",
                mean_r, mean_cpp, speedup, total_mem_mb))

    results <- rbind(results, data.frame(
      n_obs = n,
      r_time = mean_r,
      cpp_time = mean_cpp,
      speedup = speedup,
      memory_mb = total_mem_mb
    ))
  }

  cat("\n")
  return(results)
}

# Profile component functions
profile_exponential_components <- function(n = 500) {
  cat("\n=== Component-Level Profiling ===\n\n")

  # Test data
  prior_params <- c(2, 2)
  x <- matrix(rexp(50, rate = 3), ncol = 1)

  # Prior draws
  cat("Prior Draw Performance:\n")
  time_r <- system.time({
    for (i in 1:100) {
      set_use_cpp(FALSE)
      md <- ExponentialMixtureCreate(prior_params)
      PriorDraw(md, 10)
    }
  })["elapsed"]

  time_cpp <- system.time({
    for (i in 1:100) {
      exponential_prior_draw_cpp(prior_params, 10)
    }
  })["elapsed"]

  cat(sprintf("  R: %.3fs, C++: %.3fs, Speedup: %.2fx\n",
              time_r, time_cpp, time_r/time_cpp))

  # Posterior draws
  cat("\nPosterior Draw Performance:\n")
  time_r <- system.time({
    for (i in 1:100) {
      set_use_cpp(FALSE)
      md <- ExponentialMixtureCreate(prior_params)
      PosteriorDraw(md, x, 10)
    }
  })["elapsed"]

  time_cpp <- system.time({
    for (i in 1:100) {
      exponential_posterior_draw_cpp(prior_params, x, 10)
    }
  })["elapsed"]

  cat(sprintf("  R: %.3fs, C++: %.3fs, Speedup: %.2fx\n",
              time_r, time_cpp, time_r/time_cpp))
}

# Test statistical equivalence
test_statistical_equivalence <- function(n = 300, n_iter = 200) {
  cat("\n=== Statistical Equivalence Test ===\n\n")

  # Generate test data
  set.seed(999)
  data <- c(rexp(n/2, rate = 1), rexp(n/2, rate = 4))

  # Run R implementation
  set_use_cpp(FALSE)
  set.seed(123)
  dp_r <- DirichletProcessExponential(data)
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  # Run C++ implementation
  set_use_cpp(TRUE)
  set.seed(123)
  dp_cpp <- DirichletProcessExponential(data)
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  # Compare results
  cat(sprintf("Number of clusters - R: %d, C++: %d\n",
              dp_r$numberClusters, dp_cpp$numberClusters))
  cat(sprintf("Alpha - R: %.3f, C++: %.3f\n",
              dp_r$alpha, dp_cpp$alpha))

  # Check cluster sizes
  if (dp_r$numberClusters == dp_cpp$numberClusters) {
    sizes_r <- sort(table(dp_r$clusterLabels))
    sizes_cpp <- sort(table(dp_cpp$clusterLabels))

    cat("\nCluster sizes:\n")
    cat("  R:", sizes_r, "\n")
    cat("  C++:", sizes_cpp, "\n")

    if (all(sizes_r == sizes_cpp)) {
      cat("\n✓ Cluster sizes match!\n")
    } else {
      cat("\n✗ Cluster sizes differ\n")
    }
  }

  return(list(dp_r = dp_r, dp_cpp = dp_cpp))
}

# Main execution
if (interactive()) {
  cat("Running Exponential Distribution Performance Tests\n")
  cat("==============================================\n")

  # Check C++ availability
  cpp_status <- get_cpp_status()
  cat("\nC++ Implementation Status:\n")
  print(cpp_status)

  if (cpp_status$exponential) {
    # Run quick benchmark
    bench_results <- quick_exponential_benchmark()

    # Profile components
    comp_results <- profile_exponential_components(n = 500)

    # Test equivalence
    equiv_results <- test_statistical_equivalence(n = 300, n_iter = 200)

    cat("\n==============================================\n")
    cat("All tests completed successfully!\n")
  } else {
    cat("\nExponential C++ implementation not available.\n")
    cat("Please compile the package with C++ support.\n")
  }
}
