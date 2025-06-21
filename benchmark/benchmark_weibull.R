# inst/benchmarks/benchmark_weibull.R
# Comprehensive Benchmarks for Weibull Distribution
# This script performs detailed performance and memory profiling

library(dirichletprocess)
library(ggplot2)
library(dplyr)
library(tidyr)
library(microbenchmark)

# Memory tracking fallback
if (!requireNamespace("pryr", quietly = TRUE)) {
  mem_used <- function() {
    gc_info <- gc()
    (gc_info[1, 2] + gc_info[2, 2]) * 1024
  }
} else {
  mem_used <- pryr::mem_used
}

#' Generate synthetic Weibull mixture data
#'
#' @param n Number of observations
#' @param shapes Vector of shape parameters for each component
#' @param scales Vector of scale parameters for each component
#' @param weights Mixing weights (default: equal weights)
#' @return Numeric vector of observations
generate_weibull_mixture <- function(n, shapes, scales, weights = NULL) {
  k <- length(shapes)
  if (is.null(weights)) weights <- rep(1/k, k)

  # Generate cluster assignments
  clusters <- sample(1:k, n, replace = TRUE, prob = weights)

  # Generate data from each cluster
  data <- numeric(n)
  for (i in 1:k) {
    idx <- clusters == i
    data[idx] <- rweibull(sum(idx), shape = shapes[i], scale = scales[i])
  }

  return(sample(data)) # Shuffle to remove ordering
}

#' Benchmark Weibull DP across different scenarios
#'
#' @param n_obs_vec Vector of dataset sizes to test
#' @param n_iter_vec Vector of iteration counts to test
#' @param n_clusters_vec Vector of true cluster numbers to test
#' @param n_reps Number of repetitions for each scenario
#' @return Data frame with benchmark results
benchmark_weibull_comprehensive <- function(
    n_obs_vec = c(100, 500, 1000, 5000),
    n_iter_vec = c(100, 500, 1000),
    n_clusters_vec = c(2, 3, 5),
    n_reps = 5,
    priorParams = c(10, 2, 4)  # phi, alpha0, beta0
) {

  results <- list()
  scenario_id <- 1

  cat("Running comprehensive Weibull benchmarks...\n")
  cat("======================================\n")

  for (n_obs in n_obs_vec) {
    for (n_iter in n_iter_vec) {
      for (n_clusters in n_clusters_vec) {

        cat(sprintf("Scenario %d: n=%d, iter=%d, clusters=%d\n",
                    scenario_id, n_obs, n_iter, n_clusters))

        # Generate data with distinct clusters
        shapes <- seq(0.5, 3, length.out = n_clusters)
        scales <- seq(0.5, 2, length.out = n_clusters)
        data <- generate_weibull_mixture(n_obs, shapes, scales)

        # Benchmark R implementation
        cat("  Running R implementation...\n")
        set_use_cpp(FALSE)

        for (rep in 1:n_reps) {
          gc() # Clean memory

          time_r <- system.time({
            dp_r <- DirichletProcessWeibull(data, priorParams)
            dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
          })

          results[[length(results) + 1]] <- data.frame(
            scenario_id = scenario_id,
            n_obs = n_obs,
            n_iter = n_iter,
            n_true_clusters = n_clusters,
            rep = rep,
            implementation = "R",
            time = time_r["elapsed"],
            memory_mb = as.numeric(object.size(dp_r)) / 1024^2,
            n_clusters_found = dp_r$numberClusters,
            final_alpha = dp_r$alpha,
            stringsAsFactors = FALSE
          )
        }

        # Benchmark C++ implementation
        cat("  Running C++ implementation...\n")
        set_use_cpp(TRUE)

        for (rep in 1:n_reps) {
          gc() # Clean memory
          clear_memory_tracking()

          time_cpp <- system.time({
            dp_cpp <- DirichletProcessWeibull(data, priorParams)
            dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
          })

          mem_track <- get_memory_tracking()

          results[[length(results) + 1]] <- data.frame(
            scenario_id = scenario_id,
            n_obs = n_obs,
            n_iter = n_iter,
            n_true_clusters = n_clusters,
            rep = rep,
            implementation = "C++",
            time = time_cpp["elapsed"],
            memory_mb = as.numeric(object.size(dp_cpp)) / 1024^2,
            n_clusters_found = dp_cpp$numberClusters,
            final_alpha = dp_cpp$alpha,
            stringsAsFactors = FALSE
          )
        }

        scenario_id <- scenario_id + 1
      }
    }
  }

  return(bind_rows(results))
}

#' Profile memory usage for Weibull DP
#'
#' @param n_obs Number of observations
#' @param n_iter Number of MCMC iterations
#' @param priorParams Prior parameters
#' @return List with memory profiling results
profile_weibull_memory <- function(n_obs = 1000, n_iter = 200,
                                   priorParams = c(10, 2, 4)) {

  cat("\nProfiling memory usage for Weibull DP...\n")

  # Generate test data
  data <- generate_weibull_mixture(n_obs,
                                   shapes = c(1, 2, 3),
                                   scales = c(0.5, 1, 1.5))

  # Profile R implementation
  set_use_cpp(FALSE)
  mem_trace_r <- list()
  mem_trace_r$start <- mem_used()

  dp_r <- DirichletProcessWeibull(data, priorParams)
  mem_trace_r$after_init <- mem_used()

  # Track memory during fit
  mem_during_fit_r <- numeric(10)
  for (i in 1:10) {
    dp_r <- Fit(dp_r, n_iter/10, progressBar = FALSE)
    mem_during_fit_r[i] <- as.numeric(mem_used())
  }

  mem_trace_r$during_fit <- mem_during_fit_r
  mem_trace_r$final <- mem_used()

  # Profile C++ implementation
  set_use_cpp(TRUE)
  clear_memory_tracking()

  mem_trace_cpp <- list()
  mem_trace_cpp$start <- mem_used()

  dp_cpp <- DirichletProcessWeibull(data, priorParams)
  mem_trace_cpp$after_init <- mem_used()

  # Track memory during fit
  mem_during_fit_cpp <- numeric(10)
  cpp_internal_mem <- list()

  for (i in 1:10) {
    clear_memory_tracking()
    dp_cpp <- Fit(dp_cpp, n_iter/10, progressBar = FALSE)
    mem_during_fit_cpp[i] <- as.numeric(mem_used())
    cpp_internal_mem[[i]] <- get_memory_tracking()
  }

  mem_trace_cpp$during_fit <- mem_during_fit_cpp
  mem_trace_cpp$final <- mem_used()
  mem_trace_cpp$internal_tracking <- cpp_internal_mem

  return(list(
    r_memory = mem_trace_r,
    cpp_memory = mem_trace_cpp,
    summary = data.frame(
      implementation = c("R", "C++"),
      peak_memory_mb = c(
        max(unlist(mem_trace_r)) / 1024^2,
        max(unlist(mem_trace_cpp)) / 1024^2
      ),
      final_memory_mb = c(
        as.numeric(mem_trace_r$final - mem_trace_r$start) / 1024^2,
        as.numeric(mem_trace_cpp$final - mem_trace_cpp$start) / 1024^2
      )
    )
  ))
}

#' Benchmark specific algorithmic components for Weibull
#'
#' @param n_obs Number of observations
#' @param n_clusters Number of clusters to initialize
#' @param priorParams Prior parameters
#' @return Data frame with component-level timings
benchmark_weibull_components <- function(n_obs = 1000, n_clusters = 3,
                                         priorParams = c(10, 2, 4)) {

  cat("Benchmarking individual components for Weibull...\n")

  # Generate test data
  shapes <- seq(0.8, 2.5, length.out = n_clusters)
  scales <- seq(0.5, 1.5, length.out = n_clusters)
  data_vec <- generate_weibull_mixture(n_obs, shapes, scales)
  data <- matrix(data_vec, ncol = 1)

  # Initialize DP objects
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessWeibull(data_vec, priorParams)
  dp_r <- Initialise(dp_r, numInitialClusters = n_clusters)

  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessWeibull(data_vec, priorParams)
  dp_cpp <- Initialise(dp_cpp, numInitialClusters = n_clusters)

  # Benchmark individual components
  n_reps <- 100

  # Component: Likelihood calculation
  cat("  - Likelihood calculation...\n")
  time_lik_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      for (i in 1:nrow(data)) {
        Likelihood(dp_r$mixingDistribution,
                   data[i, , drop = FALSE],
                   dp_r$clusterParameters)
      }
    },
    times = n_reps
  )

  time_lik_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      for (i in 1:nrow(data)) {
        Likelihood(dp_cpp$mixingDistribution,
                   data[i, , drop = FALSE],
                   dp_cpp$clusterParameters)
      }
    },
    times = n_reps
  )

  # Component: Cluster assignment update (non-conjugate)
  cat("  - Cluster assignment update...\n")
  time_assign_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      ClusterComponentUpdate(dp_r)
    },
    times = n_reps
  )

  time_assign_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      ClusterComponentUpdate(dp_cpp)
    },
    times = n_reps
  )

  # Component: Parameter update (MH sampling)
  cat("  - Parameter update (MH)...\n")
  time_param_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      ClusterParameterUpdate(dp_r)
    },
    times = n_reps
  )

  time_param_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      ClusterParameterUpdate(dp_cpp)
    },
    times = n_reps
  )

  # Combine results
  results <- rbind(
    data.frame(
      component = "Likelihood",
      implementation = c("R", "C++"),
      mean_time_ms = c(mean(time_lik_r$time) / 1e6,
                       mean(time_lik_cpp$time) / 1e6),
      sd_time_ms = c(sd(time_lik_r$time) / 1e6,
                     sd(time_lik_cpp$time) / 1e6)
    ),
    data.frame(
      component = "ClusterAssignment",
      implementation = c("R", "C++"),
      mean_time_ms = c(mean(time_assign_r$time) / 1e6,
                       mean(time_assign_cpp$time) / 1e6),
      sd_time_ms = c(sd(time_assign_r$time) / 1e6,
                     sd(time_assign_cpp$time) / 1e6)
    ),
    data.frame(
      component = "ParameterUpdate",
      implementation = c("R", "C++"),
      mean_time_ms = c(mean(time_param_r$time) / 1e6,
                       mean(time_param_cpp$time) / 1e6),
      sd_time_ms = c(sd(time_param_r$time) / 1e6,
                     sd(time_param_cpp$time) / 1e6)
    )
  )

  return(results)
}

# Run benchmarks and generate report
run_and_report <- function() {
  # 1. Run comprehensive benchmarks
  bench_results <- benchmark_weibull_comprehensive(
    n_obs_vec = c(100, 500, 1000, 2500),
    n_iter_vec = c(100, 500),
    n_clusters_vec = c(2, 3, 5),
    n_reps = 3
  )

  # 2. Memory profiling
  memory_results <- profile_weibull_memory(n_obs = 1000, n_iter = 200)

  # 3. Component benchmarks
  component_results <- benchmark_weibull_components(n_obs = 1000, n_clusters = 3)

  # Generate summary statistics
  summary_stats <- bench_results %>%
    group_by(implementation) %>%
    summarise(
      mean_time = mean(time),
      median_time = median(time),
      mean_memory = mean(memory_mb),
      .groups = "drop"
    )

  # Calculate speedup
  speedup <- bench_results %>%
    select(scenario_id, implementation, time, rep) %>%
    pivot_wider(names_from = implementation, values_from = time) %>%
    summarise(
      mean_speedup = mean(R / `C++`, na.rm = TRUE),
      median_speedup = median(R / `C++`, na.rm = TRUE),
      min_speedup = min(R / `C++`, na.rm = TRUE),
      max_speedup = max(R / `C++`, na.rm = TRUE)
    )

  cat("\n===== WEIBULL DP BENCHMARK SUMMARY =====\n\n")

  cat("1. OVERALL PERFORMANCE:\n")
  print(summary_stats)

  cat("\n2. SPEEDUP ANALYSIS:\n")
  print(speedup)

  cat("\n3. MEMORY USAGE:\n")
  print(memory_results$summary)

  cat("\n4. COMPONENT-LEVEL PERFORMANCE:\n")
  print(component_results)

  cat("\n5. COMPONENT SPEEDUP:\n")
  comp_speedup <- component_results %>%
    select(component, implementation, mean_time_ms) %>%
    pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
    mutate(speedup = R / `C++`)
  print(comp_speedup)

  # Create visualization
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Speedup by scenario
    speedup_plot <- bench_results %>%
      select(n_obs, n_iter, implementation, time) %>%
      group_by(n_obs, n_iter, implementation) %>%
      summarise(mean_time = mean(time), .groups = "drop") %>%
      pivot_wider(names_from = implementation, values_from = mean_time) %>%
      mutate(speedup = R / `C++`) %>%
      ggplot(aes(x = n_obs, y = speedup, color = as.factor(n_iter))) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      scale_x_log10() +
      labs(
        title = "Weibull DP: C++ Speedup vs Dataset Size",
        x = "Number of Observations",
        y = "Speedup Factor (R time / C++ time)",
        color = "Iterations"
      ) +
      theme_minimal()

    print(speedup_plot)
  }

  return(list(
    benchmarks = bench_results,
    memory = memory_results,
    components = component_results,
    summary_stats = summary_stats,
    speedup = speedup
  ))
}

# Example usage:
results <- run_and_report()
