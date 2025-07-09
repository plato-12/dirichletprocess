# inst/benchmarks/benchmark_exponential.R
# Comprehensive Benchmarks for Exponential Distribution
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

#' Generate synthetic exponential mixture data
#'
#' @param n Number of observations
#' @param rates Vector of rate parameters for each component
#' @param weights Mixing weights (default: equal weights)
#' @return Numeric vector of observations
generate_exponential_mixture <- function(n, rates, weights = NULL) {
  k <- length(rates)
  if (is.null(weights)) weights <- rep(1/k, k)

  # Generate cluster assignments
  clusters <- sample(1:k, n, replace = TRUE, prob = weights)

  # Generate data from each cluster
  data <- numeric(n)
  for (i in 1:k) {
    idx <- clusters == i
    data[idx] <- rexp(sum(idx), rate = rates[i])
  }

  return(sample(data)) # Shuffle to remove ordering
}

#' Benchmark exponential DP across different scenarios
#'
#' @param n_obs_vec Vector of dataset sizes to test
#' @param n_iter_vec Vector of iteration counts to test
#' @param n_clusters_vec Vector of true cluster numbers to test
#' @param n_reps Number of repetitions for each scenario
#' @return Data frame with benchmark results
benchmark_exponential_comprehensive <- function(
    n_obs_vec = c(100, 500, 1000, 5000),
    n_iter_vec = c(100, 500, 1000),
    n_clusters_vec = c(2, 3, 5),
    n_reps = 5
) {

  results <- list()
  scenario_id <- 1

  cat("Running comprehensive exponential benchmarks...\n")
  cat("======================================\n")

  for (n_obs in n_obs_vec) {
    for (n_iter in n_iter_vec) {
      for (n_clusters in n_clusters_vec) {

        cat(sprintf("Scenario %d: n=%d, iter=%d, clusters=%d\n",
                    scenario_id, n_obs, n_iter, n_clusters))

        # Generate data with exponential mixture
        rates <- seq(0.5, 5, length.out = n_clusters)
        data <- generate_exponential_mixture(n_obs, rates)

        # Storage for results
        times_r <- numeric(n_reps)
        times_cpp <- numeric(n_reps)
        memory_r <- numeric(n_reps)
        memory_cpp <- numeric(n_reps)
        clusters_found_r <- numeric(n_reps)
        clusters_found_cpp <- numeric(n_reps)

        for (rep in 1:n_reps) {
          # R implementation
          set_use_cpp(FALSE)
          gc()
          mem_before <- as.numeric(mem_used())

          time_r <- system.time({
            dp_r <- DirichletProcessExponential(data)
            dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
          })

          mem_after <- as.numeric(mem_used())
          times_r[rep] <- time_r["elapsed"]
          memory_r[rep] <- mem_after - mem_before
          clusters_found_r[rep] <- dp_r$numberClusters

          # C++ implementation
          set_use_cpp(TRUE)
          gc()
          mem_before <- as.numeric(mem_used())

          time_cpp <- system.time({
            dp_cpp <- DirichletProcessExponential(data)
            dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
          })

          mem_after <- as.numeric(mem_used())
          times_cpp[rep] <- time_cpp["elapsed"]
          memory_cpp[rep] <- mem_after - mem_before
          clusters_found_cpp[rep] <- dp_cpp$numberClusters
        }

        # Store results
        results[[scenario_id]] <- data.frame(
          scenario_id = scenario_id,
          n_obs = n_obs,
          n_iter = n_iter,
          n_true_clusters = n_clusters,
          implementation = rep(c("R", "C++"), each = n_reps),
          time = c(times_r, times_cpp),
          memory_mb = c(memory_r, memory_cpp) / 1024^2,
          clusters_found = c(clusters_found_r, clusters_found_cpp),
          rep = rep(1:n_reps, 2)
        )

        scenario_id <- scenario_id + 1
      }
    }
  }

  cat("======================================\n")
  cat("Benchmarking complete!\n")

  return(bind_rows(results))
}

#' Detailed memory profiling for exponential DP
#'
#' @param n_obs Number of observations
#' @param n_iter Number of iterations
#' @return List with memory profiling results
profile_exponential_memory <- function(n_obs = 1000, n_iter = 500) {

  cat("Running detailed memory profiling...\n")

  # Generate test data
  data <- generate_exponential_mixture(n_obs, rates = c(1, 3, 5))

  # Profile R implementation
  cat("Profiling R implementation...\n")
  set_use_cpp(FALSE)

  # Track memory at key points
  mem_trace_r <- list()
  gc()
  mem_trace_r$start <- pryr::mem_used()

  dp_r <- DirichletProcessExponential(data)
  mem_trace_r$after_init <- pryr::mem_used()

  # Sample memory usage during fitting
  mem_during_fit_r <- numeric(10)
  iter_per_sample <- n_iter / 10

  for (i in 1:10) {
    dp_r <- Fit(dp_r, iter_per_sample, progressBar = FALSE)
    mem_during_fit_r[i] <- as.numeric(pryr::mem_used())
  }

  mem_trace_r$during_fit <- mem_during_fit_r
  mem_trace_r$final <- pryr::mem_used()

  # Profile C++ implementation
  cat("Profiling C++ implementation...\n")
  set_use_cpp(TRUE)
  clear_memory_tracking()

  mem_trace_cpp <- list()
  gc()
  mem_trace_cpp$start <- pryr::mem_used()

  dp_cpp <- DirichletProcessExponential(data)
  mem_trace_cpp$after_init <- pryr::mem_used()

  # Sample memory usage during fitting
  mem_during_fit_cpp <- numeric(10)
  cpp_internal_mem <- list()

  for (i in 1:10) {
    dp_cpp <- Fit(dp_cpp, iter_per_sample, progressBar = FALSE)
    mem_during_fit_cpp[i] <- as.numeric(pryr::mem_used())
    cpp_internal_mem[[i]] <- get_memory_tracking()
  }

  mem_trace_cpp$during_fit <- mem_during_fit_cpp
  mem_trace_cpp$final <- pryr::mem_used()
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

#' Benchmark specific algorithmic components
#'
#' @param n_obs Number of observations
#' @param n_clusters Number of clusters to initialize
#' @return Data frame with component-level timings
benchmark_exponential_components <- function(n_obs = 1000, n_clusters = 3) {

  cat("Benchmarking individual components...\n")

  # Generate test data
  data_vec <- generate_exponential_mixture(n_obs, rates = seq(0.5, 5, length.out = n_clusters))
  data <- matrix(data_vec, ncol = 1)

  # Initialize DP objects
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessExponential(data)
  dp_r <- Initialise(dp_r, numInitialClusters = n_clusters)

  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessExponential(data)
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

  # Component: Cluster assignment update
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

  # Component: Parameter update
  cat("  - Parameter update...\n")
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

#' Test statistical equivalence between R and C++ implementations
#'
#' @param n Number of observations
#' @param n_iter Number of iterations
#' @param n_chains Number of independent chains to run
#' @return List with comparison results
test_statistical_equivalence <- function(n = 500, n_iter = 1000, n_chains = 3) {

  cat("Testing statistical equivalence...\n")

  # Generate test data
  set.seed(42)
  data <- c(rexp(n/3, rate = 1),
            rexp(n/3, rate = 3),
            rexp(n/3, rate = 5))

  # Run multiple chains for each implementation
  r_results <- list()
  cpp_results <- list()

  for (chain in 1:n_chains) {
    cat(sprintf("  Chain %d/%d\n", chain, n_chains))

    # R implementation
    set_use_cpp(FALSE)
    set.seed(1000 + chain)
    dp_r <- DirichletProcessExponential(data)
    dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

    r_results[[chain]] <- list(
      n_clusters = dp_r$numberClusters,
      alpha = dp_r$alpha,
      cluster_sizes = table(dp_r$clusterLabels)
    )

    # C++ implementation
    set_use_cpp(TRUE)
    set.seed(1000 + chain)
    dp_cpp <- DirichletProcessExponential(data)
    dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

    cpp_results[[chain]] <- list(
      n_clusters = dp_cpp$numberClusters,
      alpha = dp_cpp$alpha,
      cluster_sizes = table(dp_cpp$clusterLabels)
    )
  }

  # Compare results
  r_clusters <- sapply(r_results, function(x) x$n_clusters)
  cpp_clusters <- sapply(cpp_results, function(x) x$n_clusters)

  r_alpha <- sapply(r_results, function(x) x$alpha)
  cpp_alpha <- sapply(cpp_results, function(x) x$alpha)

  cat("\nStatistical Equivalence Results:\n")
  cat("================================\n")
  cat(sprintf("Mean clusters - R: %.2f, C++: %.2f\n",
              mean(r_clusters), mean(cpp_clusters)))
  cat(sprintf("Mean alpha - R: %.3f, C++: %.3f\n",
              mean(r_alpha), mean(cpp_alpha)))

  return(list(
    r_results = r_results,
    cpp_results = cpp_results,
    comparison = data.frame(
      metric = c("n_clusters", "alpha"),
      r_mean = c(mean(r_clusters), mean(r_alpha)),
      cpp_mean = c(mean(cpp_clusters), mean(cpp_alpha)),
      r_sd = c(sd(r_clusters), sd(r_alpha)),
      cpp_sd = c(sd(cpp_clusters), sd(cpp_alpha))
    )
  ))
}

#' Profile component-level performance
#'
#' @param n Number of observations
#' @return List with profiling results
profile_exponential_components <- function(n = 1000) {

  cat("\nProfiling exponential components...\n")

  # Generate data
  data <- generate_exponential_mixture(n, rates = c(1, 3, 5))

  # Test prior and posterior draws
  prior_params <- c(2, 2)  # alpha0, beta0
  x_sample <- matrix(rexp(20, rate = 2), ncol = 1)

  # Time prior draws
  cat("  Testing prior draws...\n")
  prior_bench <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      md <- ExponentialMixtureCreate(prior_params)
      PriorDraw(md, 100)
    },
    Cpp = {
      exponential_prior_draw_cpp(prior_params, 100)
    },
    times = 100
  )

  # Time posterior draws
  cat("  Testing posterior draws...\n")
  post_bench <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      md <- ExponentialMixtureCreate(prior_params)
      PosteriorDraw(md, x_sample, 100)
    },
    Cpp = {
      exponential_posterior_draw_cpp(prior_params, x_sample, 100)
    },
    times = 100
  )

  # Time likelihood calculations
  cat("  Testing likelihood calculations...\n")
  x_vec <- seq(0.1, 5, length.out = 100)
  lambda <- 2.5

  lik_bench <- microbenchmark(
    R = dexp(x_vec, lambda),
    Cpp = exponential_likelihood_cpp(x_vec, lambda),
    times = 100
  )

  # Create summary
  summary_df <- rbind(
    data.frame(
      operation = "PriorDraw",
      R_ms = mean(prior_bench$time[prior_bench$expr == "R"]) / 1e6,
      Cpp_ms = mean(prior_bench$time[prior_bench$expr == "Cpp"]) / 1e6
    ),
    data.frame(
      operation = "PosteriorDraw",
      R_ms = mean(post_bench$time[post_bench$expr == "R"]) / 1e6,
      Cpp_ms = mean(post_bench$time[post_bench$expr == "Cpp"]) / 1e6
    ),
    data.frame(
      operation = "Likelihood",
      R_ms = mean(lik_bench$time[lik_bench$expr == "R"]) / 1e6,
      Cpp_ms = mean(lik_bench$time[lik_bench$expr == "Cpp"]) / 1e6
    )
  )

  summary_df$speedup <- summary_df$R_ms / summary_df$Cpp_ms

  return(list(
    benchmarks = list(
      prior = prior_bench,
      posterior = post_bench,
      likelihood = lik_bench
    ),
    summary = summary_df
  ))
}

#' Run benchmarks and generate report
run_and_report <- function() {
  # 1. Run comprehensive benchmarks
  bench_results <- benchmark_exponential_comprehensive(
    n_obs_vec = c(100, 1000),
    n_iter_vec = c(100, 500),
    n_clusters_vec = c(2, 3, 5),
    n_reps = 3
  )

  # 2. Memory profiling
  memory_results <- profile_exponential_memory(n_obs = 1000, n_iter = 200)

  # 3. Component benchmarks
  component_results <- benchmark_exponential_components(n_obs = 1000, n_clusters = 3)

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

  cat("\n===== EXPONENTIAL DP BENCHMARK SUMMARY =====\n\n")

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
