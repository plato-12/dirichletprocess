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

    cat(sprintf("  R:     %.3f sec (mean of %d runs)\n", mean_r, n_reps))
    cat(sprintf("  C++:   %.3f sec (mean of %d runs)\n", mean_cpp, n_reps))
    cat(sprintf("  Speedup: %.1fx\n", speedup))
    cat(sprintf("  C++ Memory: %.2f MB\n\n", total_mem_mb))

    # Store results
    results <- rbind(results, data.frame(
      n_obs = n,
      time_r = mean_r,
      time_cpp = mean_cpp,
      speedup = speedup,
      memory_mb = total_mem_mb
    ))
  }

  return(results)
}

# Detailed component profiling
profile_exponential_components <- function(n = 500) {
  cat("\n=== Component-Level Profiling ===\n\n")

  # Generate data
  set.seed(123)
  data <- c(rexp(n/3, rate = 1), rexp(n/3, rate = 3), rexp(n/3, rate = 5))

  # Initialize both implementations
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessExponential(data)
  dp_r <- Initialise(dp_r, numInitialClusters = 3)

  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessExponential(data)
  dp_cpp <- Initialise(dp_cpp, numInitialClusters = 3)

  # Time individual components
  n_reps <- 50

  # 1. Likelihood calculations
  cat("1. Likelihood Calculations:\n")

  set_use_cpp(FALSE)
  time_lik_r <- system.time({
    for (rep in 1:n_reps) {
      for (i in 1:n) {
        Likelihood(dp_r$mixingDistribution,
                   matrix(data[i], ncol = 1),
                   dp_r$clusterParameters)
      }
    }
  })["elapsed"] / n_reps

  set_use_cpp(TRUE)
  time_lik_cpp <- system.time({
    for (rep in 1:n_reps) {
      for (i in 1:n) {
        Likelihood(dp_cpp$mixingDistribution,
                   matrix(data[i], ncol = 1),
                   dp_cpp$clusterParameters)
      }
    }
  })["elapsed"] / n_reps

  cat(sprintf("  R: %.3f sec, C++: %.3f sec, Speedup: %.1fx\n",
              time_lik_r, time_lik_cpp, time_lik_r/time_lik_cpp))

  # 2. Cluster updates
  cat("\n2. Cluster Component Updates:\n")

  set_use_cpp(FALSE)
  time_cluster_r <- system.time({
    for (rep in 1:n_reps) {
      ClusterComponentUpdate(dp_r)
    }
  })["elapsed"] / n_reps

  set_use_cpp(TRUE)
  clear_memory_tracking()
  time_cluster_cpp <- system.time({
    for (rep in 1:n_reps) {
      ClusterComponentUpdate(dp_cpp)
    }
  })["elapsed"] / n_reps
  cluster_mem <- get_memory_tracking()

  cat(sprintf("  R: %.3f sec, C++: %.3f sec, Speedup: %.1fx\n",
              time_cluster_r, time_cluster_cpp, time_cluster_r/time_cluster_cpp))

  if (nrow(cluster_mem) > 0) {
    cat("  C++ Memory allocations:\n")
    print(cluster_mem)
  }

  # 3. Parameter updates
  cat("\n3. Parameter Updates:\n")

  set_use_cpp(FALSE)
  time_param_r <- system.time({
    for (rep in 1:n_reps) {
      ClusterParameterUpdate(dp_r)
    }
  })["elapsed"] / n_reps

  set_use_cpp(TRUE)
  clear_memory_tracking()
  time_param_cpp <- system.time({
    for (rep in 1:n_reps) {
      ClusterParameterUpdate(dp_cpp)
    }
  })["elapsed"] / n_reps
  param_mem <- get_memory_tracking()

  cat(sprintf("  R: %.3f sec, C++: %.3f sec, Speedup: %.1fx\n",
              time_param_r, time_param_cpp, time_param_r/time_param_cpp))

  if (nrow(param_mem) > 0) {
    cat("  C++ Memory allocations:\n")
    print(param_mem)
  }

  return(list(
    likelihood = c(r = time_lik_r, cpp = time_lik_cpp),
    cluster_update = c(r = time_cluster_r, cpp = time_cluster_cpp),
    param_update = c(r = time_param_r, cpp = time_param_cpp)
  ))
}

# Test statistical equivalence
test_statistical_equivalence <- function(n = 200, n_iter = 100) {
  cat("\n=== Testing Statistical Equivalence ===\n\n")

  # Generate test data
  set.seed(456)
  true_rates <- c(1, 3, 5)
  data <- c(
    rexp(n/3, rate = true_rates[1]),
    rexp(n/3, rate = true_rates[2]),
    rexp(n/3, rate = true_rates[3])
  )

  # Fit with R implementation
  set.seed(789)
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessExponential(data)
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  # Fit with C++ implementation
  set.seed(789)
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessExponential(data)
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  # Compare results
  cat("Number of clusters found:\n")
  cat(sprintf("  R: %d, C++: %d\n", dp_r$numberClusters, dp_cpp$numberClusters))

  # Compare cluster parameters (rates)
  cat("\nCluster rates:\n")
  cat("  R:", sort(as.numeric(dp_r$clusterParameters[[1]])), "\n")
  cat("  C++:", sort(as.numeric(dp_cpp$clusterParameters[[1]])), "\n")

  # Compare final likelihood
  cat("\nFinal likelihood:\n")
  cat(sprintf("  R: %.3f, C++: %.3f\n",
              tail(dp_r$likelihoodChain, 1),
              tail(dp_cpp$likelihoodChain, 1)))

  # Test if cluster assignments are similar
  if (dp_r$numberClusters == dp_cpp$numberClusters) {
    # Simple check: compare sorted cluster sizes
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
