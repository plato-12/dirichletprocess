# benchmark/comprehensive_performance_tests.R

library(microbenchmark)
library(ggplot2)

# Performance test configuration
PERFORMANCE_CONFIG <- list(
  sample_sizes = c(100, 500, 1000, 5000, 10000),
  iterations = c(10, 50, 100, 500),
  distributions = c("normal", "exponential", "beta", "weibull", "mvnormal", "mvnormal2"),
  n_replications = 10
)

run_comprehensive_performance_tests <- function() {
  results <- list()

  for (dist in PERFORMANCE_CONFIG$distributions) {
    cat("\nTesting", dist, "distribution...\n")

    for (n in PERFORMANCE_CONFIG$sample_sizes) {
      cat("  Sample size:", n, "\n")

      # Generate test data
      test_data <- generate_test_data(dist, n)

      for (its in PERFORMANCE_CONFIG$iterations) {
        cat("    Iterations:", its, "")

        # Benchmark R vs C++
        bench_result <- microbenchmark(
          R = {
            set_use_cpp(FALSE)
            dp <- create_dp_object(dist, test_data)
            Fit(dp, its = its)
          },
          Cpp = {
            set_use_cpp(TRUE)
            dp <- create_dp_object(dist, test_data)
            Fit(dp, its = its)
          },
          times = PERFORMANCE_CONFIG$n_replications
        )

        # Store results
        results[[length(results) + 1]] <- list(
          distribution = dist,
          sample_size = n,
          iterations = its,
          benchmark = bench_result,
          speedup = median(bench_result$time[bench_result$expr == "R"]) /
            median(bench_result$time[bench_result$expr == "Cpp"])
        )

        cat(" Speedup:", round(results[[length(results)]]$speedup, 2), "x\n")
      }
    }
  }

  # Create performance report
  create_performance_report(results)

  return(results)
}

# Scaling analysis
test_scaling_behavior <- function() {
  cat("\nTesting scaling behavior...\n")

  scaling_results <- list()

  for (dist in c("normal", "mvnormal")) {
    sample_sizes <- seq(100, 10000, by = 500)

    r_times <- numeric(length(sample_sizes))
    cpp_times <- numeric(length(sample_sizes))

    for (i in seq_along(sample_sizes)) {
      n <- sample_sizes[i]
      test_data <- generate_test_data(dist, n)

      # R timing
      set_use_cpp(FALSE)
      r_time <- system.time({
        dp <- create_dp_object(dist, test_data)
        Fit(dp, its = 50)
      })[3]
      r_times[i] <- r_time

      # C++ timing
      set_use_cpp(TRUE)
      cpp_time <- system.time({
        dp <- create_dp_object(dist, test_data)
        Fit(dp, its = 50)
      })[3]
      cpp_times[i] <- cpp_time
    }

    scaling_results[[dist]] <- data.frame(
      n = sample_sizes,
      r_time = r_times,
      cpp_time = cpp_times,
      speedup = r_times / cpp_times
    )
  }

  # Plot scaling behavior
  plot_scaling_results(scaling_results)

  return(scaling_results)
}

# Memory profiling
profile_memory_usage <- function() {
  if (!requireNamespace("profmem", quietly = TRUE)) {
    warning("profmem package not installed, skipping memory profiling")
    return(NULL)
  }

  library(profmem)

  memory_results <- list()

  for (dist in c("normal", "mvnormal")) {
    test_data <- generate_test_data(dist, n = 1000)

    # Profile R implementation
    set_use_cpp(FALSE)
    r_mem <- profmem({
      dp <- create_dp_object(dist, test_data)
      Fit(dp, its = 100)
    })

    # Profile C++ implementation
    set_use_cpp(TRUE)
    cpp_mem <- profmem({
      dp <- create_dp_object(dist, test_data)
      Fit(dp, its = 100)
    })

    memory_results[[dist]] <- list(
      r_memory = total(r_mem),
      cpp_memory = total(cpp_mem),
      memory_ratio = total(r_mem) / total(cpp_mem)
    )
  }

  return(memory_results)
}

# Quick benchmark for development
quick_benchmark <- function(dist = "normal", n = 1000, its = 100) {
  cat("\nQuick benchmark for", dist, "with n =", n, "and", its, "iterations\n")

  test_data <- generate_test_data(dist, n)

  bench_result <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      dp <- create_dp_object(dist, test_data)
      Fit(dp, its = its)
    },
    Cpp = {
      set_use_cpp(TRUE)
      dp <- create_dp_object(dist, test_data)
      Fit(dp, its = its)
    },
    times = 5
  )

  print(bench_result)

  speedup <- median(bench_result$time[bench_result$expr == "R"]) /
    median(bench_result$time[bench_result$expr == "Cpp"])

  cat("\nSpeedup factor:", round(speedup, 2), "x\n")

  return(bench_result)
}

# Performance test by algorithm type
test_algorithm_performance <- function() {
  cat("\nTesting performance by algorithm type...\n")

  # Conjugate distributions (Algorithm 4)
  conjugate_dists <- c("normal", "exponential")

  # Non-conjugate distributions (Algorithm 8)
  nonconjugate_dists <- c("beta", "weibull")

  results <- list()

  for (dist_type in c("conjugate", "nonconjugate")) {
    dists <- if(dist_type == "conjugate") conjugate_dists else nonconjugate_dists

    for (dist in dists) {
      test_data <- generate_test_data(dist, 500)

      bench <- microbenchmark(
        R = {
          set_use_cpp(FALSE)
          dp <- create_dp_object(dist, test_data)
          Fit(dp, its = 100)
        },
        Cpp = {
          set_use_cpp(TRUE)
          dp <- create_dp_object(dist, test_data)
          Fit(dp, its = 100)
        },
        times = 10
      )

      results[[dist]] <- list(
        type = dist_type,
        benchmark = bench,
        speedup = median(bench$time[bench$expr == "R"]) /
          median(bench$time[bench$expr == "Cpp"])
      )
    }
  }

  # Compare conjugate vs non-conjugate speedups
  conjugate_speedups <- sapply(results[conjugate_dists], `[[`, "speedup")
  nonconjugate_speedups <- sapply(results[nonconjugate_dists], `[[`, "speedup")

  cat("\nConjugate average speedup:", round(mean(conjugate_speedups), 2), "x\n")
  cat("Non-conjugate average speedup:", round(mean(nonconjugate_speedups), 2), "x\n")

  return(results)
}
