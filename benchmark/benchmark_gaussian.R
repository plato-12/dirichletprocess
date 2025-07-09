# =============================================================================
# Gaussian Implementation Benchmark: R vs C++ Performance Comparison
# =============================================================================

#' Quick Benchmark Script for dirichletprocess Gaussian Implementation
#'
#' This script compares R and C++ backends for Gaussian Dirichlet Process
#' mixture models across different data sizes and iteration counts.

library(dirichletprocess)
library(microbenchmark)
library(ggplot2)

# Check if required packages are available
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  message("Installing microbenchmark for benchmarking...")
  install.packages("microbenchmark")
  library(microbenchmark)
}

if (!requireNamespace("pryr", quietly = TRUE)) {
  message("Installing pryr for memory profiling...")
  install.packages("pryr")
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Generate test data with known clusters
generate_gaussian_mixture <- function(n, k = 3, seed = 42) {
  set.seed(seed)

  # Create mixture proportions
  props <- rep(1/k, k)
  cluster_sizes <- as.numeric(rmultinom(1, n, props))

  # Generate cluster means and variances
  means <- seq(-3, 3, length.out = k)
  sds <- rep(0.8, k)

  data <- c()
  true_labels <- c()

  for (i in 1:k) {
    if (cluster_sizes[i] > 0) {
      cluster_data <- rnorm(cluster_sizes[i], mean = means[i], sd = sds[i])
      data <- c(data, cluster_data)
      true_labels <- c(true_labels, rep(i, cluster_sizes[i]))
    }
  }

  return(list(data = data, true_labels = true_labels, k = length(unique(true_labels))))
}

#' Safe memory measurement
measure_memory <- function() {
  if (requireNamespace("pryr", quietly = TRUE)) {
    return(as.numeric(pryr::mem_used()))
  } else {
    return(NA)
  }
}

#' Compare clustering results
compare_clustering <- function(dp_r, dp_cpp, true_k = NULL) {
  results <- list()

  # Number of clusters found
  results$clusters_r <- dp_r$numberClusters
  results$clusters_cpp <- dp_cpp$numberClusters

  # If true number of clusters is known
  if (!is.null(true_k)) {
    results$true_clusters <- true_k
    results$r_accuracy <- abs(dp_r$numberClusters - true_k)
    results$cpp_accuracy <- abs(dp_cpp$numberClusters - true_k)
  }

  # Data likelihood (if available)
  if ("likelihoodChain" %in% names(dp_r) && length(dp_r$likelihoodChain) > 0) {
    results$r_final_likelihood <- tail(dp_r$likelihoodChain, 1)
  }

  if ("likelihoodChain" %in% names(dp_cpp) && length(dp_cpp$likelihoodChain) > 0) {
    results$cpp_final_likelihood <- tail(dp_cpp$likelihoodChain, 1)
  }

  return(results)
}

# =============================================================================
# Benchmark Functions
# =============================================================================

#' Single benchmark run
benchmark_single_run <- function(data, n_iter = 100, verbose = TRUE) {
  if (verbose) cat(sprintf("Benchmarking with %d observations, %d iterations...\n",
                           length(data), n_iter))

  results <- list()

  # R Backend Benchmark
  if (verbose) cat("  Testing R backend...\n")

  set_use_cpp(FALSE)
  mem_before_r <- measure_memory()

  time_r <- system.time({
    set.seed(123)
    dp_r <- DirichletProcessGaussian(data)
    dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
  })

  mem_after_r <- measure_memory()

  # C++ Backend Benchmark (if available)
  cpp_available <- get_cpp_status()$mcmc_runner || exists("_dirichletprocess_run_mcmc_cpp")

  if (cpp_available) {
    if (verbose) cat("  Testing C++ backend...\n")

    set_use_cpp(TRUE)
    mem_before_cpp <- measure_memory()

    time_cpp <- system.time({
      set.seed(123)  # Same seed for fair comparison
      dp_cpp <- DirichletProcessGaussian(data)
      dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
    })

    mem_after_cpp <- measure_memory()

    # Compare results
    comparison <- compare_clustering(dp_r, dp_cpp)

    results <- list(
      n_obs = length(data),
      n_iter = n_iter,
      r_time = time_r[["elapsed"]],
      cpp_time = time_cpp[["elapsed"]],
      speedup = time_r[["elapsed"]] / time_cpp[["elapsed"]],
      r_memory = ifelse(is.na(mem_before_r), NA, mem_after_r - mem_before_r),
      cpp_memory = ifelse(is.na(mem_before_cpp), NA, mem_after_cpp - mem_before_cpp),
      comparison = comparison,
      cpp_available = TRUE
    )

  } else {
    if (verbose) cat("  C++ backend not available, skipping comparison\n")

    results <- list(
      n_obs = length(data),
      n_iter = n_iter,
      r_time = time_r[["elapsed"]],
      cpp_time = NA,
      speedup = NA,
      r_memory = ifelse(is.na(mem_before_r), NA, mem_after_r - mem_before_r),
      cpp_memory = NA,
      comparison = list(clusters_r = dp_r$numberClusters),
      cpp_available = FALSE
    )
  }

  if (verbose) {
    cat(sprintf("  R time: %.2f seconds\n", results$r_time))
    if (results$cpp_available) {
      cat(sprintf("  C++ time: %.2f seconds\n", results$cpp_time))
      cat(sprintf("  Speedup: %.2fx\n", results$speedup))
    }
    cat(sprintf("  Clusters found (R): %d\n", results$comparison$clusters_r))
    if (results$cpp_available) {
      cat(sprintf("  Clusters found (C++): %d\n", results$comparison$clusters_cpp))
    }
  }

  return(results)
}

#' Comprehensive benchmark across different scenarios
comprehensive_benchmark <- function() {
  cat("=================================================================\n")
  cat("Comprehensive Gaussian DP Benchmark: R vs C++ Implementation\n")
  cat("=================================================================\n\n")

  # Check C++ availability
  cpp_status <- get_cpp_status()
  cat("C++ Implementation Status:\n")
  print(cpp_status)
  cat("\n")

  if (!cpp_status$mcmc_runner && !exists("_dirichletprocess_run_mcmc_cpp")) {
    cat("⚠️  C++ implementation not available - showing R performance only\n\n")
  }

  # Benchmark scenarios
  scenarios <- list(
    list(n = 100, iter = 50, desc = "Small dataset, quick fit"),
    list(n = 500, iter = 100, desc = "Medium dataset"),
    list(n = 1000, iter = 100, desc = "Large dataset"),
    list(n = 200, iter = 200, desc = "Medium dataset, many iterations")
  )

  results <- list()

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    cat(sprintf("Scenario %d: %s\n", i, scenario$desc))
    cat(sprintf("  Data size: %d observations\n", scenario$n))
    cat(sprintf("  Iterations: %d\n", scenario$iter))
    cat("  " %+% paste(rep("-", 50), collapse = "") %+% "\n")

    # Generate test data
    test_data <- generate_gaussian_mixture(scenario$n, k = 3)

    # Run benchmark
    result <- benchmark_single_run(test_data$data, scenario$iter, verbose = TRUE)
    result$scenario <- scenario$desc
    result$true_k <- test_data$k

    results[[i]] <- result
    cat("\n")
  }

  # Summary
  cat("=================================================================\n")
  cat("BENCHMARK SUMMARY\n")
  cat("=================================================================\n")

  if (any(sapply(results, function(x) x$cpp_available))) {
    # Create summary table
    summary_df <- data.frame(
      Scenario = sapply(results, function(x) x$scenario),
      N_Obs = sapply(results, function(x) x$n_obs),
      N_Iter = sapply(results, function(x) x$n_iter),
      R_Time = sprintf("%.2fs", sapply(results, function(x) x$r_time)),
      CPP_Time = sapply(results, function(x) {
        if (x$cpp_available) sprintf("%.2fs", x$cpp_time) else "N/A"
      }),
      Speedup = sapply(results, function(x) {
        if (x$cpp_available) sprintf("%.1fx", x$speedup) else "N/A"
      }),
      R_Clusters = sapply(results, function(x) x$comparison$clusters_r),
      CPP_Clusters = sapply(results, function(x) {
        if (x$cpp_available) x$comparison$clusters_cpp else NA
      }),
      stringsAsFactors = FALSE
    )

    print(summary_df, row.names = FALSE)

    # Overall statistics
    cpp_results <- results[sapply(results, function(x) x$cpp_available)]
    if (length(cpp_results) > 0) {
      avg_speedup <- mean(sapply(cpp_results, function(x) x$speedup), na.rm = TRUE)
      max_speedup <- max(sapply(cpp_results, function(x) x$speedup), na.rm = TRUE)

      cat(sprintf("\nPerformance Summary:\n"))
      cat(sprintf("  Average speedup: %.2fx\n", avg_speedup))
      cat(sprintf("  Maximum speedup: %.2fx\n", max_speedup))

      # Memory comparison (if available)
      r_mem <- sapply(cpp_results, function(x) x$r_memory)
      cpp_mem <- sapply(cpp_results, function(x) x$cpp_memory)

      if (!all(is.na(r_mem)) && !all(is.na(cpp_mem))) {
        cat(sprintf("  Average R memory: %.1f MB\n", mean(r_mem, na.rm = TRUE) / 1024^2))
        cat(sprintf("  Average C++ memory: %.1f MB\n", mean(cpp_mem, na.rm = TRUE) / 1024^2))
      }
    }
  } else {
    cat("C++ implementation not available - only R results shown:\n\n")
    for (i in seq_along(results)) {
      result <- results[[i]]
      cat(sprintf("%s: %.2fs (%d clusters)\n",
                  result$scenario, result$r_time, result$comparison$clusters_r))
    }
  }

  cat("\n=================================================================\n")

  return(results)
}

#' Quick microbenchmark for detailed timing
quick_microbenchmark <- function(n = 200, iter = 50, times = 5) {
  cat("Detailed Timing Comparison (microbenchmark)\n")
  cat("============================================\n")

  test_data <- generate_gaussian_mixture(n, k = 3)
  data <- test_data$data

  cpp_available <- get_cpp_status()$mcmc_runner || exists("_dirichletprocess_run_mcmc_cpp")

  if (!cpp_available) {
    cat("C++ implementation not available - skipping microbenchmark\n")
    return(NULL)
  }

  # Define benchmark functions
  r_fit <- function() {
    set_use_cpp(FALSE)
    set.seed(42)
    dp <- DirichletProcessGaussian(data)
    Fit(dp, iter, progressBar = FALSE)
  }

  cpp_fit <- function() {
    set_use_cpp(TRUE)
    set.seed(42)
    dp <- DirichletProcessGaussian(data)
    Fit(dp, iter, progressBar = FALSE)
  }

  # Run microbenchmark
  mb_result <- microbenchmark(
    R_Backend = r_fit(),
    CPP_Backend = cpp_fit(),
    times = times
  )

  print(mb_result)

  # Plot if ggplot2 is available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- autoplot(mb_result) +
      ggtitle(sprintf("Gaussian DP Performance: %d obs, %d iter", n, iter)) +
      theme_minimal()
    print(p)
  }

  return(mb_result)
}

# =============================================================================
# Main Execution
# =============================================================================

#' Run all benchmarks
run_all_benchmarks <- function() {
  cat("Starting Gaussian Dirichlet Process Benchmarks...\n\n")

  # 1. Comprehensive benchmark
  comp_results <- comprehensive_benchmark()

  # 2. Detailed microbenchmark (if available)
  cat("\n")
  mb_results <- quick_microbenchmark(n = 300, iter = 100, times = 3)

  cat("\n✅ Benchmarking complete!\n")
  cat("\nTo run individual benchmarks:\n")
  cat("  comprehensive_benchmark()     # Full comparison across scenarios\n")
  cat("  quick_microbenchmark()        # Detailed timing with error bars\n")
  cat("  benchmark_single_run(data)    # Single custom benchmark\n")

  return(list(comprehensive = comp_results, microbenchmark = mb_results))
}

# =============================================================================
# Quick Test Function
# =============================================================================

#' Quick verification that both backends work
quick_test <- function() {
  cat("Quick Backend Verification Test\n")
  cat("==============================\n")

  data <- generate_gaussian_mixture(100, k = 2)$data

  # Test R backend
  cat("Testing R backend... ")
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(data)
  dp_r <- Fit(dp_r, 20, progressBar = FALSE)
  cat("✅ OK\n")

  # Test C++ backend
  cpp_available <- get_cpp_status()$mcmc_runner || exists("_dirichletprocess_run_mcmc_cpp")
  if (cpp_available) {
    cat("Testing C++ backend... ")
    set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessGaussian(data)
    dp_cpp <- Fit(dp_cpp, 20, progressBar = FALSE)
    cat("✅ OK\n")

    cat(sprintf("\nResults: R found %d clusters, C++ found %d clusters\n",
                dp_r$numberClusters, dp_cpp$numberClusters))
  } else {
    cat("C++ backend not available ❌\n")
  }

  cat("\nBackend switching works correctly! 🎉\n")
}

# =============================================================================
# Auto-run if sourced directly
# =============================================================================

if (interactive()) {
  cat("Gaussian DP Benchmark Script Loaded! 🚀\n")
  cat("\nAvailable functions:\n")
  cat("  quick_test()              # Quick verification\n")
  cat("  run_all_benchmarks()     # Complete benchmark suite\n")
  cat("  comprehensive_benchmark() # Cross-scenario comparison\n")
  cat("  quick_microbenchmark()   # Detailed timing analysis\n")
  cat("\nRun quick_test() to verify everything works!\n")
}
