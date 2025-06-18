# Benchmark Beta Distribution Implementation
# Compare R vs C++ performance, memory usage, and scaling

library(dirichletprocess)
library(microbenchmark)
library(ggplot2)
library(dplyr)
library(tidyr)

# Source additional benchmarking utilities if available
if (file.exists("inst/benchmarks/utils.R")) {
  source("inst/benchmarks/utils.R")
}

#' Generate Beta mixture data
#'
#' @param n Number of observations
#' @param shapes List of (alpha, beta) pairs for each component
#' @param weights Mixture weights (default: equal)
#' @param seed Random seed
#' @return Vector of observations
generate_beta_mixture <- function(n, shapes = list(c(2, 8), c(8, 2)),
                                  weights = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  k <- length(shapes)
  if (is.null(weights)) weights <- rep(1/k, k)

  # Generate component assignments
  components <- sample(1:k, n, replace = TRUE, prob = weights)

  # Generate data
  data <- numeric(n)
  for (i in 1:n) {
    comp <- components[i]
    data[i] <- rbeta(1, shapes[[comp]][1], shapes[[comp]][2])
  }

  return(data)
}

#' Run comprehensive Beta benchmark
#'
#' @param n_obs_vec Vector of observation counts to test
#' @param n_iter Number of MCMC iterations
#' @param n_reps Number of replications per scenario
#' @param shapes List of shape parameters for data generation
#' @return Data frame with benchmark results
benchmark_beta_comprehensive <- function(
    n_obs_vec = c(50, 100, 200, 500, 1000),
    n_iter = 100,
    n_reps = 5,
    shapes = list(
      single = list(c(5, 5)),
      two_cluster = list(c(2, 8), c(8, 2)),
      three_cluster = list(c(2, 8), c(5, 5), c(8, 2))
    )) {

  cat("======================================\n")
  cat("Beta Distribution Comprehensive Benchmark\n")
  cat("======================================\n")
  cat("Testing configurations:\n")
  cat("- Observations:", paste(n_obs_vec, collapse = ", "), "\n")
  cat("- Iterations:", n_iter, "\n")
  cat("- Replications:", n_reps, "\n")
  cat("- Scenarios:", paste(names(shapes), collapse = ", "), "\n")
  cat("======================================\n\n")

  results <- list()
  scenario_id <- 1

  for (scenario_name in names(shapes)) {
    cat("\nScenario:", scenario_name, "\n")

    for (n_obs in n_obs_vec) {
      cat(sprintf("  n = %d: ", n_obs))

      # Generate data
      data <- generate_beta_mixture(n_obs, shapes[[scenario_name]], seed = 123)

      # Storage for this configuration
      times_r <- numeric(n_reps)
      times_cpp <- numeric(n_reps)
      memory_r <- numeric(n_reps)
      memory_cpp <- numeric(n_reps)
      clusters_found_r <- numeric(n_reps)
      clusters_found_cpp <- numeric(n_reps)
      final_alpha_r <- numeric(n_reps)
      final_alpha_cpp <- numeric(n_reps)

      # Run replications
      for (rep in 1:n_reps) {
        cat(".")

        # R implementation
        set_use_cpp(FALSE)
        gc()
        mem_before <- pryr::mem_used()

        time_r <- system.time({
          dp_r <- DirichletProcessBeta(data, verbose = FALSE)
          dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
        })["elapsed"]

        mem_after <- pryr::mem_used()

        times_r[rep] <- time_r
        memory_r[rep] <- as.numeric(mem_after - mem_before)
        clusters_found_r[rep] <- dp_r$numberClusters
        final_alpha_r[rep] <- dp_r$alpha

        # C++ implementation (if available)
        if (can_use_cpp(DirichletProcessBeta(data[1:10]))) {
          set_use_cpp(TRUE)
          gc()
          mem_before <- pryr::mem_used()

          time_cpp <- system.time({
            dp_cpp <- DirichletProcessBeta(data, verbose = FALSE)
            dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
          })["elapsed"]

          mem_after <- pryr::mem_used()

          times_cpp[rep] <- time_cpp
          memory_cpp[rep] <- as.numeric(mem_after - mem_before)
          clusters_found_cpp[rep] <- dp_cpp$numberClusters
          final_alpha_cpp[rep] <- dp_cpp$alpha
        } else {
          times_cpp[rep] <- NA
          memory_cpp[rep] <- NA
          clusters_found_cpp[rep] <- NA
          final_alpha_cpp[rep] <- NA
        }
      }

      cat(" Done\n")

      # Store results
      results[[scenario_id]] <- data.frame(
        scenario = scenario_name,
        n_obs = n_obs,
        n_iter = n_iter,
        n_true_clusters = length(shapes[[scenario_name]]),
        implementation = rep(c("R", "C++"), each = n_reps),
        time = c(times_r, times_cpp),
        memory_mb = c(memory_r, memory_cpp) / 1024^2,
        clusters_found = c(clusters_found_r, clusters_found_cpp),
        final_alpha = c(final_alpha_r, final_alpha_cpp),
        rep = rep(1:n_reps, 2)
      )

      scenario_id <- scenario_id + 1
    }
  }

  cat("======================================\n")
  cat("Benchmarking complete!\n")

  return(bind_rows(results))
}

#' Profile Beta memory usage in detail
#'
#' @param n_obs Number of observations
#' @param n_iter Number of iterations
#' @return List with memory profiling results
profile_beta_memory <- function(n_obs = 1000, n_iter = 500) {

  cat("Running detailed memory profiling...\n")

  # Generate test data - three component mixture
  data <- generate_beta_mixture(n_obs,
                                shapes = list(c(2, 8), c(5, 5), c(8, 2)),
                                weights = c(0.3, 0.4, 0.3))

  # Profile R implementation
  cat("Profiling R implementation...\n")
  set_use_cpp(FALSE)

  mem_trace_r <- list()
  gc()
  mem_trace_r$start <- pryr::mem_used()

  dp_r <- DirichletProcessBeta(data, verbose = FALSE)
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

  # Profile C++ implementation if available
  if (can_use_cpp(DirichletProcessBeta(data[1:10]))) {
    cat("Profiling C++ implementation...\n")
    set_use_cpp(TRUE)

    if (exists("clear_memory_tracking")) {
      clear_memory_tracking()
    }

    mem_trace_cpp <- list()
    gc()
    mem_trace_cpp$start <- pryr::mem_used()

    dp_cpp <- DirichletProcessBeta(data, verbose = FALSE)
    mem_trace_cpp$after_init <- pryr::mem_used()

    mem_during_fit_cpp <- numeric(10)
    cpp_internal_mem <- list()

    for (i in 1:10) {
      dp_cpp <- Fit(dp_cpp, iter_per_sample, progressBar = FALSE)
      mem_during_fit_cpp[i] <- as.numeric(pryr::mem_used())

      if (exists("get_memory_tracking")) {
        cpp_internal_mem[[i]] <- get_memory_tracking()
      }
    }

    mem_trace_cpp$during_fit <- mem_during_fit_cpp
    mem_trace_cpp$final <- pryr::mem_used()
    mem_trace_cpp$internal_tracking <- cpp_internal_mem
  } else {
    mem_trace_cpp <- NULL
  }

  return(list(
    r_memory = mem_trace_r,
    cpp_memory = mem_trace_cpp,
    summary = data.frame(
      implementation = c("R", if (!is.null(mem_trace_cpp)) "C++" else NULL),
      peak_memory_mb = c(
        max(unlist(mem_trace_r)) / 1024^2,
        if (!is.null(mem_trace_cpp)) max(unlist(mem_trace_cpp[names(mem_trace_cpp) != "internal_tracking"])) / 1024^2 else NULL
      ),
      final_memory_mb = c(
        as.numeric(mem_trace_r$final - mem_trace_r$start) / 1024^2,
        if (!is.null(mem_trace_cpp)) as.numeric(mem_trace_cpp$final - mem_trace_cpp$start) / 1024^2 else NULL
      )
    )
  ))
}

#' Benchmark Beta algorithmic components
#'
#' @param n_obs Number of observations
#' @param n_clusters Number of clusters to initialize
#' @return Data frame with component-level timings
benchmark_beta_components <- function(n_obs = 1000, n_clusters = 3) {

  cat("Benchmarking individual components...\n")

  # Generate test data
  data_vec <- generate_beta_mixture(n_obs,
                                    shapes = lapply(1:n_clusters, function(i) c(i*2, 10-i*2)))
  data <- matrix(data_vec, ncol = 1)

  # Initialize DP objects
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessBeta(data, verbose = FALSE)
  dp_r <- Initialise(dp_r, numInitialClusters = n_clusters)

  results <- list()

  if (can_use_cpp(dp_r)) {
    set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessBeta(data, verbose = FALSE)
    dp_cpp <- Initialise(dp_cpp, numInitialClusters = n_clusters)

    # Benchmark individual components
    n_reps <- 100

    # Component: Likelihood calculation
    cat("  - Likelihood calculation...\n")

    # Single point likelihood
    single_point_results <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        Likelihood(dp_r$mixingDistribution,
                   data[1, , drop = FALSE],
                   dp_r$clusterParameters)
      },
      Cpp = {
        set_use_cpp(TRUE)
        Likelihood(dp_cpp$mixingDistribution,
                   data[1, , drop = FALSE],
                   dp_cpp$clusterParameters)
      },
      times = n_reps
    )

    # All points likelihood
    all_points_results <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        for (i in 1:nrow(data)) {
          Likelihood(dp_r$mixingDistribution,
                     data[i, , drop = FALSE],
                     dp_r$clusterParameters)
        }
      },
      Cpp = {
        set_use_cpp(TRUE)
        for (i in 1:nrow(data)) {
          Likelihood(dp_cpp$mixingDistribution,
                     data[i, , drop = FALSE],
                     dp_cpp$clusterParameters)
        }
      },
      times = 10
    )

    # Component: Cluster label updates
    cat("  - Cluster label updates...\n")

    cluster_update_results <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        ClusterComponentUpdate(dp_r)
      },
      Cpp = {
        set_use_cpp(TRUE)
        ClusterComponentUpdate(dp_cpp)
      },
      times = n_reps
    )

    # Component: Parameter updates
    cat("  - Parameter updates...\n")

    param_update_results <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        ClusterParameterUpdate(dp_r)
      },
      Cpp = {
        set_use_cpp(TRUE)
        ClusterParameterUpdate(dp_cpp)
      },
      times = n_reps
    )

    # Component: Prior draws
    cat("  - Prior draws...\n")

    prior_draw_results <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        PriorDraw(dp_r$mixingDistribution, 10)
      },
      Cpp = {
        set_use_cpp(TRUE)
        PriorDraw(dp_cpp$mixingDistribution, 10)
      },
      times = n_reps
    )

    # Component: Posterior draws (MH steps)
    cat("  - Posterior draws (MH)...\n")

    cluster_data <- data[dp_r$clusterLabels == 1, , drop = FALSE]

    posterior_draw_results <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        PosteriorDraw(dp_r$mixingDistribution, cluster_data, 1)
      },
      Cpp = {
        set_use_cpp(TRUE)
        PosteriorDraw(dp_cpp$mixingDistribution, cluster_data, 1)
      },
      times = 20  # Fewer reps as this is expensive
    )

    # Compile results
    results <- list(
      single_likelihood = single_point_results,
      all_likelihood = all_points_results,
      cluster_update = cluster_update_results,
      param_update = param_update_results,
      prior_draw = prior_draw_results,
      posterior_draw = posterior_draw_results
    )
  }

  return(results)
}

#' Test Beta statistical equivalence between R and C++
#'
#' @param n Number of observations
#' @param n_iter Number of MCMC iterations
#' @return List with comparison results
test_beta_statistical_equivalence <- function(n = 200, n_iter = 500) {

  cat("\n=== Testing Statistical Equivalence ===\n\n")

  # Generate test data from known mixture
  set.seed(456)
  true_shapes <- list(c(2, 8), c(5, 5), c(8, 2))
  true_weights <- c(0.3, 0.4, 0.3)
  data <- generate_beta_mixture(n, true_shapes, true_weights)

  # Fit with R implementation
  cat("Fitting with R implementation...\n")
  set.seed(789)
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessBeta(data, verbose = FALSE)
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  results <- list(R = dp_r)

  # Fit with C++ implementation if available
  if (can_use_cpp(DirichletProcessBeta(data[1:10]))) {
    cat("Fitting with C++ implementation...\n")
    set.seed(789)
    set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessBeta(data, verbose = FALSE)
    dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

    results$Cpp <- dp_cpp

    # Compare results
    cat("\nComparison Results:\n")
    cat("==================\n")

    cat(sprintf("Number of clusters - R: %d, C++: %d\n",
                dp_r$numberClusters, dp_cpp$numberClusters))

    cat(sprintf("Final alpha - R: %.3f, C++: %.3f\n",
                dp_r$alpha, dp_cpp$alpha))

    # Compare cluster parameters
    cat("\nCluster Parameters:\n")
    cat("R implementation:\n")
    for (k in 1:dp_r$numberClusters) {
      cat(sprintf("  Cluster %d: mu=%.3f, tau=%.3f (n=%d)\n",
                  k, dp_r$clusterParameters[[1]][k],
                  dp_r$clusterParameters[[2]][k],
                  dp_r$pointsPerCluster[k]))
    }

    cat("\nC++ implementation:\n")
    for (k in 1:dp_cpp$numberClusters) {
      cat(sprintf("  Cluster %d: mu=%.3f, tau=%.3f (n=%d)\n",
                  k, dp_cpp$clusterParameters[[1]][k],
                  dp_cpp$clusterParameters[[2]][k],
                  dp_cpp$pointsPerCluster[k]))
    }

    # Statistical tests for equivalence
    if (length(dp_r$likelihoodChain) > 0 && length(dp_cpp$likelihoodChain) > 0) {
      # Compare final likelihoods
      final_lik_r <- tail(dp_r$likelihoodChain, 100)
      final_lik_cpp <- tail(dp_cpp$likelihoodChain, 100)

      ks_test <- ks.test(final_lik_r, final_lik_cpp)
      cat(sprintf("\nKS test for likelihood chains: p-value = %.3f\n",
                  ks_test$p.value))
    }
  }

  return(results)
}

#' Create Beta benchmark visualizations
#'
#' @param results Data frame from benchmark_beta_comprehensive
#' @param output_dir Directory to save plots
create_beta_benchmark_plots <- function(results, output_dir = "inst/benchmarks/plots") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Remove NA values
  results_clean <- results %>%
    filter(!is.na(time))

  # 1. Timing comparison plot
  p1 <- ggplot(results_clean, aes(x = n_obs, y = time, color = implementation)) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
    facet_wrap(~ scenario, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Beta DP Performance: R vs C++",
         x = "Number of observations",
         y = "Time (seconds)",
         color = "Implementation") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "beta_timing_comparison.png"), p1,
         width = 10, height = 6, dpi = 300)

  # 2. Speedup plot
  speedup_data <- results_clean %>%
    group_by(scenario, n_obs, rep) %>%
    summarise(speedup = time[implementation == "R"] / time[implementation == "C++"],
              .groups = "drop") %>%
    filter(!is.na(speedup), is.finite(speedup))

  p2 <- ggplot(speedup_data, aes(x = n_obs, y = speedup)) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~ scenario) +
    scale_x_log10() +
    labs(title = "C++ Speedup over R Implementation",
         x = "Number of observations",
         y = "Speedup factor",
         caption = "Values > 1 indicate C++ is faster") +
    theme_minimal()

  ggsave(file.path(output_dir, "beta_speedup.png"), p2,
         width = 10, height = 6, dpi = 300)

  # 3. Memory usage plot
  p3 <- ggplot(results_clean, aes(x = n_obs, y = memory_mb, color = implementation)) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
    facet_wrap(~ scenario) +
    scale_x_log10() +
    labs(title = "Memory Usage: R vs C++",
         x = "Number of observations",
         y = "Memory (MB)",
         color = "Implementation") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "beta_memory_usage.png"), p3,
         width = 10, height = 6, dpi = 300)

  # 4. Clustering accuracy plot
  p4 <- ggplot(results_clean, aes(x = implementation, y = clusters_found)) +
    geom_boxplot(aes(fill = implementation)) +
    geom_hline(aes(yintercept = n_true_clusters), linetype = "dashed", color = "red") +
    facet_grid(scenario ~ n_obs) +
    labs(title = "Clustering Performance",
         x = "Implementation",
         y = "Number of clusters found",
         caption = "Red line indicates true number of clusters") +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "beta_clustering_accuracy.png"), p4,
         width = 12, height = 8, dpi = 300)

  cat("\nPlots saved to:", output_dir, "\n")
}

#' Generate Beta benchmark report
#'
#' @param results Benchmark results
#' @param memory_profile Memory profiling results
#' @param component_results Component benchmark results
#' @param output_file Path for output report
generate_beta_benchmark_report <- function(results, memory_profile,
                                           component_results,
                                           output_file = "beta_benchmark_report.txt") {

  sink(output_file)

  cat("=====================================\n")
  cat("Beta Distribution Benchmark Report\n")
  cat("=====================================\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # Summary statistics
  summary_stats <- results %>%
    filter(!is.na(time)) %>%
    group_by(implementation, scenario, n_obs) %>%
    summarise(
      mean_time = mean(time),
      sd_time = sd(time),
      mean_memory = mean(memory_mb),
      mean_clusters = mean(clusters_found),
      .groups = "drop"
    )

  cat("Performance Summary\n")
  cat("==================\n\n")

  # Calculate overall speedup
  speedup_summary <- results %>%
    filter(!is.na(time)) %>%
    group_by(scenario, n_obs, rep) %>%
    summarise(
      speedup = time[implementation == "R"] / time[implementation == "C++"],
      .groups = "drop"
    ) %>%
    filter(!is.na(speedup), is.finite(speedup)) %>%
    group_by(scenario) %>%
    summarise(
      mean_speedup = mean(speedup),
      min_speedup = min(speedup),
      max_speedup = max(speedup),
      .groups = "drop"
    )

  cat("C++ Speedup by Scenario:\n")
  print(speedup_summary)
  cat("\n")

  # Memory usage summary
  cat("\nMemory Usage Summary:\n")
  cat("====================\n")
  if (!is.null(memory_profile$summary)) {
    print(memory_profile$summary)
  }

  # Component-level performance
  cat("\n\nComponent-Level Performance:\n")
  cat("===========================\n")

  if (length(component_results) > 0) {
    for (comp_name in names(component_results)) {
      if (!is.null(component_results[[comp_name]])) {
        cat("\n", comp_name, ":\n", sep = "")
        comp_summary <- summary(component_results[[comp_name]])
        print(comp_summary[, c("expr", "mean", "median")])
      }
    }
  }

  # Detailed timing table
  cat("\n\nDetailed Timing Results:\n")
  cat("=======================\n")

  timing_table <- results %>%
    filter(!is.na(time)) %>%
    group_by(implementation, scenario, n_obs) %>%
    summarise(
      mean_time = mean(time),
      sd_time = sd(time),
      min_time = min(time),
      max_time = max(time),
      .groups = "drop"
    ) %>%
    arrange(scenario, n_obs, implementation)

  print(as.data.frame(timing_table), row.names = FALSE)

  sink()
  cat("\nReport saved to:", output_file, "\n")
}

# Main execution function
if (sys.nframe() == 0) {
  cat("Running Beta Distribution Benchmarks\n")
  cat("===================================\n\n")

  # Check if C++ is available
  cpp_available <- can_use_cpp(DirichletProcessBeta(rbeta(10, 2, 2)))
  cat("C++ implementation available:", cpp_available, "\n\n")

  # Run comprehensive benchmark
  results <- benchmark_beta_comprehensive(
    n_obs_vec = c(50, 100, 200, 500, 1000),
    n_iter = 100,
    n_reps = 5
  )

  # Memory profiling
  memory_profile <- profile_beta_memory(n_obs = 500, n_iter = 200)

  # Component benchmarks
  component_results <- benchmark_beta_components(n_obs = 500)

  # Statistical equivalence
  equivalence_results <- test_beta_statistical_equivalence(n = 200, n_iter = 200)

  # Generate plots
  create_beta_benchmark_plots(results)

  # Generate report
  generate_beta_benchmark_report(results, memory_profile,
                                 component_results)

  # Print summary
  cat("\n\nBenchmark Summary\n")
  cat("=================\n")

  if (cpp_available) {
    avg_speedup <- results %>%
      filter(!is.na(time)) %>%
      group_by(n_obs, rep) %>%
      summarise(speedup = time[implementation == "R"] / time[implementation == "C++"],
                .groups = "drop") %>%
      filter(!is.na(speedup), is.finite(speedup)) %>%
      pull(speedup) %>%
      mean()

    cat(sprintf("Average C++ speedup: %.1fx\n", avg_speedup))
  }

  cat("\nBenchmarking complete!\n")
}
