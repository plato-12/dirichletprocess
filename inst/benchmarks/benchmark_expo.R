#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive Performance Benchmark for dirichletprocess Package
# Comparing R vs C++ Implementation for Exponential Mixture Models
# =============================================================================

# Load required libraries
library(dirichletprocess)
library(microbenchmark)
library(profmem)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

# Set random seed for reproducibility
set.seed(42)

# =============================================================================
# Configuration Parameters
# =============================================================================

# Data generation parameters
TRUE_CLUSTERS <- 3
TRUE_RATES <- c(1.0, 3.0, 7.0)  # Exponential rate parameters
CLUSTER_PROPORTIONS <- c(0.3, 0.5, 0.2)

# Benchmark parameters
DATA_SIZES <- c(100, 500, 1000, 2500, 5000)
N_ITERATIONS <- c(100, 500, 1000)  # MCMC iterations to test
N_BENCHMARK_RUNS <- 10  # Number of times to repeat each benchmark
WARMUP_RUNS <- 2  # Warmup runs before timing

# Output settings
OUTPUT_DIR <- "benchmark_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# =============================================================================
# Helper Functions
# =============================================================================

#' Generate synthetic exponential mixture data
#' @param n Number of data points
#' @param rates Vector of exponential rate parameters
#' @param proportions Mixing proportions (must sum to 1)
#' @return Numeric vector of generated data
generate_exponential_mixture <- function(n, rates, proportions) {
  k <- length(rates)
  cluster_assignments <- sample(1:k, n, replace = TRUE, prob = proportions)

  data <- numeric(n)
  for (i in 1:k) {
    idx <- cluster_assignments == i
    data[idx] <- rexp(sum(idx), rate = rates[i])
  }

  return(data)
}

#' Run complete MCMC sampler
#' @param dp DirichletProcess object
#' @param n_iter Number of iterations
#' @param use_cpp Whether to use C++ implementation
#' @return Updated DirichletProcess object
run_mcmc_complete <- function(dp, n_iter, use_cpp = FALSE) {
  # Set implementation ONCE at the beginning
  dirichletprocess:::set_use_cpp(use_cpp)

  # Run the complete MCMC chain
  for (i in 1:n_iter) {
    dp <- dirichletprocess:::ClusterComponentUpdate(dp)
    dp <- dirichletprocess:::ClusterParameterUpdate(dp)
    dp <- dirichletprocess:::UpdateAlpha(dp)
  }

  return(dp)
}

#' Benchmark execution time with proper isolation
#' @param data_sizes Vector of data sizes to test
#' @param n_iter Number of MCMC iterations
#' @param n_runs Number of benchmark runs
#' @return Data frame with benchmark results
benchmark_execution_time <- function(data_sizes, n_iter, n_runs) {
  results <- list()

  cat("\n=== Execution Time Benchmarking ===\n")

  for (size in data_sizes) {
    cat(sprintf("\nBenchmarking with n=%d, iterations=%d...\n", size, n_iter))

    # Generate data
    data <- generate_exponential_mixture(size, TRUE_RATES, CLUSTER_PROPORTIONS)

    # Create base DP object
    dp_base <- DirichletProcessExponential(data)

    # Benchmark R implementation
    cat("  Testing R implementation...\n")
    dirichletprocess:::set_use_cpp(FALSE)  # Ensure R mode

    # Warmup for R
    for (i in 1:WARMUP_RUNS) {
      dp_temp <- dp_base
      run_mcmc_complete(dp_temp, 10, use_cpp = FALSE)
    }

    bench_r <- microbenchmark(
      R = {
        dp_temp <- dp_base
        run_mcmc_complete(dp_temp, n_iter, use_cpp = FALSE)
      },
      times = n_runs
    )

    # Benchmark C++ implementation
    cat("  Testing C++ implementation...\n")
    dirichletprocess:::set_use_cpp(TRUE)  # Ensure C++ mode

    # Warmup for C++
    for (i in 1:WARMUP_RUNS) {
      dp_temp <- dp_base
      run_mcmc_complete(dp_temp, 10, use_cpp = TRUE)
    }

    bench_cpp <- microbenchmark(
      Cpp = {
        dp_temp <- dp_base
        run_mcmc_complete(dp_temp, n_iter, use_cpp = TRUE)
      },
      times = n_runs
    )

    # Reset to default
    dirichletprocess:::set_use_cpp(FALSE)

    # Combine results
    results[[length(results) + 1]] <- data.frame(
      data_size = size,
      iterations = n_iter,
      implementation = c(rep("R", n_runs), rep("Cpp", n_runs)),
      time_ms = c(bench_r$time / 1e6, bench_cpp$time / 1e6)
    )

    # Print summary
    cat(sprintf("  R: mean=%.1f ms, median=%.1f ms\n",
                mean(bench_r$time) / 1e6, median(bench_r$time) / 1e6))
    cat(sprintf("  C++: mean=%.1f ms, median=%.1f ms\n",
                mean(bench_cpp$time) / 1e6, median(bench_cpp$time) / 1e6))
    cat(sprintf("  Speedup: %.2fx\n",
                mean(bench_r$time) / mean(bench_cpp$time)))
  }

  return(bind_rows(results))
}

#' Profile memory usage
#' @param data_sizes Vector of data sizes to test
#' @param n_iter Number of MCMC iterations
#' @return Data frame with memory profiling results
profile_memory_usage <- function(data_sizes, n_iter) {
  results <- list()

  cat("\n=== Memory Usage Profiling ===\n")

  for (size in data_sizes) {
    cat(sprintf("\nProfiling memory with n=%d, iterations=%d...\n", size, n_iter))

    # Generate data
    data <- generate_exponential_mixture(size, TRUE_RATES, CLUSTER_PROPORTIONS)

    # Profile R implementation
    dirichletprocess:::set_use_cpp(FALSE)
    dp_r <- DirichletProcessExponential(data)
    gc()  # Clean memory before profiling

    prof_r <- profmem({
      run_mcmc_complete(dp_r, n_iter, use_cpp = FALSE)
    })

    # Profile C++ implementation
    dirichletprocess:::set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessExponential(data)
    gc()  # Clean memory before profiling

    prof_cpp <- profmem({
      run_mcmc_complete(dp_cpp, n_iter, use_cpp = TRUE)
    })

    # Reset to default
    dirichletprocess:::set_use_cpp(FALSE)

    # Extract memory statistics
    mem_r <- sum(prof_r$bytes, na.rm = TRUE) / 1024^2  # Convert to MB
    mem_cpp <- sum(prof_cpp$bytes, na.rm = TRUE) / 1024^2

    results[[length(results) + 1]] <- data.frame(
      data_size = size,
      iterations = n_iter,
      implementation = c("R", "Cpp"),
      memory_mb = c(mem_r, mem_cpp),
      allocations = c(length(prof_r$bytes), length(prof_cpp$bytes))
    )

    cat(sprintf("  R: %.2f MB (%d allocations)\n", mem_r, length(prof_r$bytes)))
    cat(sprintf("  C++: %.2f MB (%d allocations)\n", mem_cpp, length(prof_cpp$bytes)))
  }

  return(bind_rows(results))
}

#' Direct performance comparison (similar to unit test)
#' @param data_size Size of data
#' @param n_iter Number of iterations
#' @return Data frame with timing comparison
direct_performance_comparison <- function(data_size = 1000, n_iter = 100) {
  cat(sprintf("\n=== Direct Performance Comparison (n=%d, iter=%d) ===\n",
              data_size, n_iter))

  # Generate data
  data <- generate_exponential_mixture(data_size, TRUE_RATES, CLUSTER_PROPORTIONS)

  # Time R implementation
  dirichletprocess:::set_use_cpp(FALSE)
  dp_r <- DirichletProcessExponential(data)

  time_r <- system.time({
    for (i in 1:n_iter) {
      dp_r <- dirichletprocess:::ClusterComponentUpdate(dp_r)
      dp_r <- dirichletprocess:::ClusterParameterUpdate(dp_r)
      dp_r <- dirichletprocess:::UpdateAlpha(dp_r)
    }
  })[3]  # elapsed time

  # Time C++ implementation
  dirichletprocess:::set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessExponential(data)

  time_cpp <- system.time({
    for (i in 1:n_iter) {
      dp_cpp <- dirichletprocess:::ClusterComponentUpdate(dp_cpp)
      dp_cpp <- dirichletprocess:::ClusterParameterUpdate(dp_cpp)
      dp_cpp <- dirichletprocess:::UpdateAlpha(dp_cpp)
    }
  })[3]  # elapsed time

  # Reset to default
  dirichletprocess:::set_use_cpp(FALSE)

  cat(sprintf("R implementation: %.3f seconds\n", time_r))
  cat(sprintf("C++ implementation: %.3f seconds\n", time_cpp))
  cat(sprintf("Speedup: %.2fx\n", time_r / time_cpp))

  return(data.frame(
    method = c("R", "C++"),
    time_seconds = c(time_r, time_cpp),
    speedup = c(1, time_r / time_cpp)
  ))
}

#' Verify that R and C++ implementations produce similar results
#' @param data Numeric vector of data
#' @param n_iter Number of MCMC iterations
#' @return Logical indicating if results are similar
verify_implementation_consistency <- function(data, n_iter = 100) {
  cat("\nVerifying R and C++ implementation consistency...\n")

  # Test R implementation
  set.seed(123)
  dirichletprocess:::set_use_cpp(FALSE)
  dp_r <- DirichletProcessExponential(data)
  dp_r <- run_mcmc_complete(dp_r, n_iter, use_cpp = FALSE)

  # Test C++ implementation
  set.seed(123)
  dirichletprocess:::set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessExponential(data)
  dp_cpp <- run_mcmc_complete(dp_cpp, n_iter, use_cpp = TRUE)

  # Reset
  dirichletprocess:::set_use_cpp(FALSE)

  # Compare key outputs
  n_clusters_r <- length(unique(dp_r$clusterLabels))
  n_clusters_cpp <- length(unique(dp_cpp$clusterLabels))

  # Check if number of clusters is similar (allowing for MCMC variability)
  clusters_similar <- abs(n_clusters_r - n_clusters_cpp) <= 2

  # Check if alpha values are similar
  alpha_similar <- abs(dp_r$alpha - dp_cpp$alpha) / dp_r$alpha < 0.1

  cat(sprintf("  R implementation: %d clusters, alpha = %.3f\n", n_clusters_r, dp_r$alpha))
  cat(sprintf("  C++ implementation: %d clusters, alpha = %.3f\n", n_clusters_cpp, dp_cpp$alpha))
  cat(sprintf("  Consistency check: %s\n",
              ifelse(clusters_similar && alpha_similar, "PASSED", "FAILED")))

  return(clusters_similar && alpha_similar)
}

#' Create visualization plots
#' @param time_results Data frame with timing results
#' @param memory_results Data frame with memory results
create_benchmark_plots <- function(time_results, memory_results) {
  # Theme for all plots
  theme_set(theme_minimal(base_size = 12) +
              theme(legend.position = "bottom",
                    plot.title = element_text(hjust = 0.5, face = "bold")))

  # 1. Execution time by data size
  p1 <- time_results %>%
    group_by(data_size, implementation) %>%
    summarise(
      mean_time = mean(time_ms),
      se_time = sd(time_ms) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = data_size, y = mean_time, color = implementation)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                  width = 0.05 * max(time_results$data_size)) +
    scale_x_continuous(trans = "log10", breaks = unique(time_results$data_size)) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = c("R" = "#E41A1C", "Cpp" = "#377EB8")) +
    labs(title = "Execution Time vs Data Size",
         x = "Data Size (log scale)",
         y = "Time (ms, log scale)",
         color = "Implementation") +
    annotation_logticks()

  # 2. Speedup factor
  p2 <- time_results %>%
    group_by(data_size, implementation) %>%
    summarise(mean_time = mean(time_ms), .groups = "drop") %>%
    pivot_wider(names_from = implementation, values_from = mean_time) %>%
    mutate(speedup = R / Cpp) %>%
    ggplot(aes(x = data_size, y = speedup)) +
    geom_line(size = 1.2, color = "#4DAF4A") +
    geom_point(size = 3, color = "#4DAF4A") +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(trans = "log10", breaks = unique(time_results$data_size)) +
    labs(title = "C++ Speedup Factor",
         x = "Data Size (log scale)",
         y = "Speedup (R time / C++ time)") +
    annotation_logticks(sides = "b")

  # 3. Memory usage comparison
  p3 <- memory_results %>%
    ggplot(aes(x = data_size, y = memory_mb, fill = implementation)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_x_continuous(trans = "log10", breaks = unique(memory_results$data_size)) +
    scale_fill_manual(values = c("R" = "#E41A1C", "Cpp" = "#377EB8")) +
    labs(title = "Memory Usage Comparison",
         x = "Data Size (log scale)",
         y = "Memory Usage (MB)",
         fill = "Implementation") +
    annotation_logticks(sides = "b")

  # 4. Time distribution boxplot
  p4 <- time_results %>%
    filter(data_size %in% c(min(data_size), median(unique(data_size)), max(data_size))) %>%
    ggplot(aes(x = factor(data_size), y = time_ms, fill = implementation)) +
    geom_boxplot(alpha = 0.8) +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = c("R" = "#E41A1C", "Cpp" = "#377EB8")) +
    labs(title = "Execution Time Distribution",
         x = "Data Size",
         y = "Time (ms, log scale)",
         fill = "Implementation") +
    annotation_logticks(sides = "l")

  # Combine plots
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

#' Generate benchmark report
#' @param time_results Data frame with timing results
#' @param memory_results Data frame with memory results
generate_report <- function(time_results, memory_results) {
  cat("\n=== BENCHMARK SUMMARY REPORT ===\n\n")

  # Overall speedup statistics
  speedup_stats <- time_results %>%
    group_by(data_size, implementation) %>%
    summarise(mean_time = mean(time_ms), .groups = "drop") %>%
    pivot_wider(names_from = implementation, values_from = mean_time) %>%
    mutate(speedup = R / Cpp)

  cat("C++ Speedup Factors:\n")
  print(speedup_stats %>% select(data_size, speedup) %>%
          mutate(speedup = round(speedup, 2)))

  cat(sprintf("\nAverage speedup: %.2fx\n", mean(speedup_stats$speedup)))
  cat(sprintf("Speedup range: %.2fx - %.2fx\n",
              min(speedup_stats$speedup), max(speedup_stats$speedup)))

  # Memory efficiency
  mem_efficiency <- memory_results %>%
    pivot_wider(names_from = implementation, values_from = memory_mb) %>%
    mutate(memory_reduction_pct = (R - Cpp) / R * 100)

  cat("\nMemory Usage Comparison:\n")
  print(mem_efficiency %>%
          select(data_size, iterations, memory_reduction_pct) %>%
          mutate(memory_reduction_pct = round(memory_reduction_pct, 1)))

  # Scalability analysis
  cat("\nScalability Analysis:\n")

  # Fit power law: time ~ a * n^b
  r_data <- speedup_stats %>% filter(!is.na(R))
  cpp_data <- speedup_stats %>% filter(!is.na(Cpp))

  if (nrow(r_data) > 2 && nrow(cpp_data) > 2) {
    r_model <- lm(log(R) ~ log(data_size), data = r_data)
    cpp_model <- lm(log(Cpp) ~ log(data_size), data = cpp_data)

    cat(sprintf("R implementation scaling: O(n^%.2f)\n", coef(r_model)[2]))
    cat(sprintf("C++ implementation scaling: O(n^%.2f)\n", coef(cpp_model)[2]))
  }
}

# =============================================================================
# Main Execution
# =============================================================================

main <- function() {
  cat("Starting comprehensive benchmark of dirichletprocess package\n")
  cat("============================================================\n")

  # Check if C++ implementation is available
  tryCatch({
    dirichletprocess:::set_use_cpp(TRUE)
    test_data <- rexp(10)
    test_dp <- DirichletProcessExponential(test_data)
    dirichletprocess:::set_use_cpp(FALSE)
    cat("✓ C++ implementation detected and working\n")
  }, error = function(e) {
    stop("C++ implementation not available. Please ensure the package was compiled with C++ support.")
  })

  # Run direct performance comparison first (like unit test)
  direct_comparison <- direct_performance_comparison(data_size = 1000, n_iter = 100)

  # Verify implementation consistency
  test_data <- generate_exponential_mixture(500, TRUE_RATES, CLUSTER_PROPORTIONS)
  if (!verify_implementation_consistency(test_data, n_iter = 200)) {
    warning("R and C++ implementations may not be producing consistent results!")
  }

  # Run execution time benchmarks
  time_results <- list()

  # Only run smaller test sizes if direct comparison shows expected speedup
  if (direct_comparison$speedup[2] > 1.5) {
    for (n_iter in N_ITERATIONS) {
      cat(sprintf("\n--- Running benchmarks with %d iterations ---\n", n_iter))
      results <- benchmark_execution_time(DATA_SIZES, n_iter, N_BENCHMARK_RUNS)
      time_results[[length(time_results) + 1]] <- results
    }
  } else {
    cat("\nWARNING: Direct comparison shows unexpected results. Running limited benchmarks...\n")
    results <- benchmark_execution_time(c(100, 500, 1000), 100, 5)
    time_results[[1]] <- results
  }

  time_results_df <- bind_rows(time_results)

  # Run memory profiling (use middle iteration count)
  memory_results <- profile_memory_usage(
    unique(time_results_df$data_size)[1:3],  # First 3 sizes only
    min(N_ITERATIONS)
  )

  # Save raw results
  write.csv(time_results_df,
            file.path(OUTPUT_DIR, "execution_time_results.csv"),
            row.names = FALSE)
  write.csv(memory_results,
            file.path(OUTPUT_DIR, "memory_usage_results.csv"),
            row.names = FALSE)
  write.csv(direct_comparison,
            file.path(OUTPUT_DIR, "direct_comparison_results.csv"),
            row.names = FALSE)

  # Create visualizations
  if (nrow(time_results_df) > 0) {
    pdf(file.path(OUTPUT_DIR, "benchmark_plots.pdf"), width = 12, height = 10)

    # Plot for available data
    create_benchmark_plots(time_results_df, memory_results)

    dev.off()

    # Generate summary report
    generate_report(time_results_df, memory_results)
  }

  # Save session info for reproducibility
  sink(file.path(OUTPUT_DIR, "session_info.txt"))
  cat("Benchmark completed on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  print(sessionInfo())
  sink()

  cat("\n\nBenchmark complete! Results saved to:", OUTPUT_DIR, "\n")
}

# Run the benchmark
if (!interactive()) {
  main()
} else {
  cat("Script loaded. Run main() to execute the benchmark.\n")
}
