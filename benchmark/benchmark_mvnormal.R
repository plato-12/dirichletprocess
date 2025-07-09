# benchmark/benchmark_mvnormal.R
# Comprehensive benchmarking for MVNormal Distribution
# Compares R vs C++ implementations

library(dirichletprocess)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Quick MVNormal benchmark for interactive testing
#'
#' @param n_obs Vector of observation counts to test
#' @param dimensions Vector of data dimensions to test
#' @param n_iter Number of MCMC iterations
#' @param n_reps Number of repetitions for timing
#' @export
quick_mvnormal_benchmark <- function(n_obs = c(100, 250, 500),
                                     dimensions = c(2, 3, 5),
                                     n_iter = 100,
                                     n_reps = 3) {

  cat("\n=== Quick MVNormal Distribution Benchmark ===\n\n")

  results <- data.frame()

  for (d in dimensions) {
    for (n in n_obs) {
      cat(sprintf("Testing d = %d, n = %d:\n", d, n))

      # Generate test data with 2 well-separated clusters
      set.seed(42)
      mu1 <- rep(-3, d)
      mu2 <- rep(3, d)
      Sigma <- diag(d) * 0.5

      data1 <- rmvnorm(n/2, mu1, Sigma)
      data2 <- rmvnorm(n/2, mu2, Sigma)
      y <- rbind(data1, data2)

      # Set up priors
      g0Priors <- list(
        mu0 = rep(0, d),
        Lambda = diag(d),
        kappa0 = 1,
        nu = d + 2
      )

      # Benchmark R implementation
      times_r <- numeric(n_reps)
      for (i in 1:n_reps) {
        set_use_cpp(FALSE)
        time_r <- system.time({
          dp_r <- DirichletProcessMvnormal(y, g0Priors)
          dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
        })
        times_r[i] <- time_r["elapsed"]
      }

      # Benchmark C++ implementation
      times_cpp <- numeric(n_reps)
      clusters_cpp <- numeric(n_reps)

      for (i in 1:n_reps) {
        set_use_cpp(TRUE)
        time_cpp <- system.time({
          dp_cpp <- DirichletProcessMvnormal(y, g0Priors)
          dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
        })
        times_cpp[i] <- time_cpp["elapsed"]
        clusters_cpp[i] <- dp_cpp$numberClusters
      }

      # Calculate statistics
      mean_r <- mean(times_r)
      mean_cpp <- mean(times_cpp)
      speedup <- mean_r / mean_cpp

      # Print results
      cat(sprintf("  R: %.3fs, C++: %.3fs, Speedup: %.1fx, Clusters: %.1f\n",
                  mean_r, mean_cpp, speedup, mean(clusters_cpp)))

      # Store results
      results <- rbind(results, data.frame(
        n_obs = n,
        dimension = d,
        implementation = "R",
        time = times_r,
        memory_mb = NA,
        clusters_found = NA
      ))

      results <- rbind(results, data.frame(
        n_obs = n,
        dimension = d,
        implementation = "C++",
        time = times_cpp,
        memory_mb = NA,
        clusters_found = clusters_cpp
      ))
    }
  }

  return(results)
}

#' Profile MVNormal components
#'
#' @param n Number of observations
#' @param d Data dimension
#' @param n_clusters Number of true clusters
#' @export
profile_mvnormal_components <- function(n = 500, d = 3, n_clusters = 3) {

  cat("\n=== MVNormal Component Profiling ===\n\n")

  # Generate clustered data
  set.seed(123)
  y <- matrix(0, n, d)
  n_per_cluster <- n / n_clusters

  for (k in 1:n_clusters) {
    start_idx <- floor((k-1) * n_per_cluster) + 1
    end_idx <- floor(k * n_per_cluster)

    # Spread clusters out
    mu_k <- rep(k * 5 - 10, d)
    y[start_idx:end_idx, ] <- rmvnorm(end_idx - start_idx + 1, mu_k, diag(d) * 0.5)
  }

  g0Priors <- list(
    mu0 = rep(0, d),
    Lambda = diag(d),
    kappa0 = 1,
    nu = d + 2
  )

  # Profile R implementation
  cat("Profiling R implementation...\n")
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessMvnormal(y, g0Priors)
  dp_r <- Initialise(dp_r)

  # Time individual components
  n_timing_reps <- 100

  time_cluster_r <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- ClusterComponentUpdate(dp_r)
    }
  })["elapsed"] / n_timing_reps * 1000  # Convert to ms

  time_param_r <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- ClusterParameterUpdate(dp_r)
    }
  })["elapsed"] / n_timing_reps * 1000

  # Profile C++ implementation
  cat("Profiling C++ implementation...\n")
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessMvnormal(y, g0Priors)
  dp_cpp <- Initialise(dp_cpp)

  time_cluster_cpp <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- ClusterComponentUpdate.mvnormal.cpp(dp_cpp)
    }
  })["elapsed"] / n_timing_reps * 1000

  time_param_cpp <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- ClusterParameterUpdate.mvnormal.cpp(dp_cpp)
    }
  })["elapsed"] / n_timing_reps * 1000

  # Create results dataframe
  results <- data.frame(
    component = rep(c("ClusterAssignment", "ParameterUpdate"), 2),
    implementation = rep(c("R", "C++"), each = 2),
    mean_time_ms = c(time_cluster_r, time_param_r, time_cluster_cpp, time_param_cpp)
  )

  # Print summary
  cat("\nComponent timings (milliseconds):\n")
  print(results)

  cat("\nSpeedup by component:\n")
  speedup_cluster <- time_cluster_r / time_cluster_cpp
  speedup_param <- time_param_r / time_param_cpp
  cat(sprintf("  Cluster Assignment: %.1fx\n", speedup_cluster))
  cat(sprintf("  Parameter Update: %.1fx\n", speedup_param))

  return(results)
}

#' Test statistical equivalence between implementations
#'
#' @param n Number of observations
#' @param d Data dimension
#' @param n_iter Number of MCMC iterations
#' @param seed Random seed
#' @export
test_mvnormal_statistical_equivalence <- function(n = 200, d = 2, n_iter = 500, seed = 456) {

  cat("\n=== Testing Statistical Equivalence ===\n\n")

  set.seed(seed)

  # Generate data with known structure
  mu1 <- rep(-2, d)
  mu2 <- rep(2, d)
  Sigma1 <- diag(d) * 0.5
  Sigma2 <- diag(d) * 0.8

  y1 <- rmvnorm(n/2, mu1, Sigma1)
  y2 <- rmvnorm(n/2, mu2, Sigma2)
  y <- rbind(y1, y2)

  g0Priors <- list(
    mu0 = rep(0, d),
    Lambda = diag(d) * 2,
    kappa0 = 0.5,
    nu = d + 2
  )

  # Run R implementation
  set_use_cpp(FALSE)
  set.seed(789)
  dp_r <- DirichletProcessMvnormal(y, g0Priors)
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  # Run C++ implementation
  set_use_cpp(TRUE)
  set.seed(789)
  dp_cpp <- DirichletProcessMvnormal(y, g0Priors)
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  # Compare results
  cat("Number of clusters:\n")
  cat(sprintf("  R: %d\n", dp_r$numberClusters))
  cat(sprintf("  C++: %d\n", dp_cpp$numberClusters))

  # Compare alpha values (concentration parameter)
  cat("\nConcentration parameter (alpha):\n")
  cat(sprintf("  R: %.3f\n", dp_r$alpha))
  cat(sprintf("  C++: %.3f\n", dp_cpp$alpha))

  # Compare posterior means (last 100 iterations)
  alpha_r <- tail(dp_r$alphaChain, 100)
  alpha_cpp <- tail(dp_cpp$alphaChain, 100)

  cat("\nPosterior alpha statistics:\n")
  cat(sprintf("  R: mean=%.3f, sd=%.3f\n", mean(alpha_r), sd(alpha_r)))
  cat(sprintf("  C++: mean=%.3f, sd=%.3f\n", mean(alpha_cpp), sd(alpha_cpp)))

  # Statistical test for equivalence
  if (length(alpha_r) == length(alpha_cpp)) {
    ks_test <- ks.test(alpha_r, alpha_cpp)
    cat(sprintf("\nKS test p-value: %.3f\n", ks_test$p.value))
  }

  return(list(dp_r = dp_r, dp_cpp = dp_cpp))
}

#' Comprehensive MVNormal benchmarking
#'
#' @param n_obs_vec Vector of observation counts
#' @param dim_vec Vector of dimensions
#' @param n_iter_vec Vector of iteration counts
#' @param n_clusters_vec Vector of true cluster counts
#' @param n_reps Number of repetitions
#' @export
benchmark_mvnormal_comprehensive <- function(
    n_obs_vec = c(100, 250, 500, 1000),
    dim_vec = c(2, 3, 5),
    n_iter_vec = c(100, 250),
    n_clusters_vec = c(2, 3, 5),
    n_reps = 3) {

  cat("\n=== Comprehensive MVNormal Benchmarking ===\n\n")

  total_scenarios <- length(n_obs_vec) * length(dim_vec) *
    length(n_iter_vec) * length(n_clusters_vec)
  cat(sprintf("Testing %d scenarios with %d reps each...\n", total_scenarios, n_reps))

  results <- data.frame()
  scenario_id <- 0

  for (n_obs in n_obs_vec) {
    for (d in dim_vec) {
      for (n_iter in n_iter_vec) {
        for (n_true_clusters in n_clusters_vec) {

          scenario_id <- scenario_id + 1
          cat(sprintf("\nScenario %d/%d: n=%d, d=%d, iter=%d, k=%d\n",
                      scenario_id, total_scenarios, n_obs, d, n_iter, n_true_clusters))

          # Generate data with specified clusters
          set.seed(scenario_id)
          y <- matrix(0, n_obs, d)
          n_per_cluster <- n_obs / n_true_clusters

          for (k in 1:n_true_clusters) {
            start_idx <- floor((k-1) * n_per_cluster) + 1
            end_idx <- min(floor(k * n_per_cluster), n_obs)

            # Spread clusters out more for higher dimensions
            mu_k <- rep((k - (n_true_clusters+1)/2) * 4, d)
            Sigma_k <- diag(d) * runif(1, 0.3, 0.8)

            y[start_idx:end_idx, ] <- rmvnorm(end_idx - start_idx + 1, mu_k, Sigma_k)
          }

          # Set up priors
          g0Priors <- list(
            mu0 = rep(0, d),
            Lambda = diag(d),
            kappa0 = 1,
            nu = d + 2
          )

          # Run benchmarks
          for (rep in 1:n_reps) {
            # R implementation
            set_use_cpp(FALSE)
            mem_before <- gc(reset = TRUE)

            time_r <- system.time({
              dp_r <- DirichletProcessMvnormal(y, g0Priors)
              dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
            })

            mem_after <- gc()
            memory_r <- sum(mem_after[,6] - mem_before[,6])

            results <- rbind(results, data.frame(
              scenario_id = scenario_id,
              n_obs = n_obs,
              dimension = d,
              n_iter = n_iter,
              n_true_clusters = n_true_clusters,
              rep = rep,
              implementation = "R",
              time = time_r["elapsed"],
              memory_mb = memory_r,
              clusters_found = dp_r$numberClusters
            ))

            # C++ implementation
            set_use_cpp(TRUE)
            mem_before <- gc(reset = TRUE)

            time_cpp <- system.time({
              dp_cpp <- DirichletProcessMvnormal(y, g0Priors)
              dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
            })

            mem_after <- gc()
            memory_cpp <- sum(mem_after[,6] - mem_before[,6])

            results <- rbind(results, data.frame(
              scenario_id = scenario_id,
              n_obs = n_obs,
              dimension = d,
              n_iter = n_iter,
              n_true_clusters = n_true_clusters,
              rep = rep,
              implementation = "C++",
              time = time_cpp["elapsed"],
              memory_mb = memory_cpp,
              clusters_found = dp_cpp$numberClusters
            ))
          }

          # Print progress
          current_results <- results %>%
            filter(scenario_id == !!scenario_id) %>%
            group_by(implementation) %>%
            summarise(mean_time = mean(time), .groups = "drop")

          speedup <- current_results$mean_time[current_results$implementation == "R"] /
            current_results$mean_time[current_results$implementation == "C++"]

          cat(sprintf("  Mean times - R: %.3fs, C++: %.3fs, Speedup: %.1fx\n",
                      current_results$mean_time[current_results$implementation == "R"],
                      current_results$mean_time[current_results$implementation == "C++"],
                      speedup))
        }
      }
    }
  }

  return(results)
}

#' Create visualization dashboard for MVNormal benchmarks
#'
#' @param bench_results Results from benchmark_mvnormal_comprehensive
#' @export
visualize_mvnormal_benchmarks <- function(bench_results) {

  # Set color scheme
  impl_colors <- c("R" = "#E41A1C", "C++" = "#377EB8")

  # 1. Speedup by dimension and data size
  speedup_data <- bench_results %>%
    group_by(n_obs, dimension, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop") %>%
    pivot_wider(names_from = implementation, values_from = mean_time) %>%
    mutate(speedup = R / `C++`)

  p_speedup <- ggplot(speedup_data, aes(x = n_obs, y = speedup, color = factor(dimension))) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "C++ Speedup by Data Size and Dimension",
      x = "Number of Observations",
      y = "Speedup Factor (R time / C++ time)",
      color = "Dimension"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5)

  # 2. Scaling behavior
  scaling_data <- bench_results %>%
    filter(n_iter == 100) %>%
    group_by(n_obs, dimension, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop")

  p_scaling <- ggplot(scaling_data, aes(x = n_obs, y = mean_time,
                                        color = implementation,
                                        linetype = factor(dimension))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Computation Time Scaling",
      x = "Number of Observations",
      y = "Time (seconds)",
      color = "Implementation",
      linetype = "Dimension"
    ) +
    theme_minimal()

  # 3. Dimension impact
  dim_impact <- bench_results %>%
    filter(n_obs == 500) %>%
    group_by(dimension, implementation) %>%
    summarise(
      mean_time = mean(time),
      se_time = sd(time) / sqrt(n()),
      .groups = "drop"
    )

  p_dimension <- ggplot(dim_impact, aes(x = dimension, y = mean_time,
                                        fill = implementation)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                  position = position_dodge(0.8), width = 0.25) +
    scale_fill_manual(values = impl_colors) +
    labs(
      title = "Impact of Data Dimension (n=500)",
      x = "Dimension",
      y = "Time (seconds)",
      fill = "Implementation"
    ) +
    theme_minimal()

  # 4. Clustering accuracy
  cluster_accuracy <- bench_results %>%
    group_by(n_true_clusters, dimension, implementation) %>%
    summarise(
      mean_found = mean(clusters_found),
      sd_found = sd(clusters_found),
      .groups = "drop"
    )

  p_clusters <- ggplot(cluster_accuracy,
                       aes(x = n_true_clusters, y = mean_found,
                           color = implementation, shape = factor(dimension))) +
    geom_point(size = 3, position = position_dodge(0.3)) +
    geom_errorbar(aes(ymin = mean_found - sd_found, ymax = mean_found + sd_found),
                  position = position_dodge(0.3), width = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Clustering Accuracy",
      x = "True Number of Clusters",
      y = "Mean Clusters Found",
      color = "Implementation",
      shape = "Dimension"
    ) +
    theme_minimal()

  # Combine plots
  dashboard <- (p_speedup | p_scaling) / (p_dimension | p_clusters) +
    plot_annotation(
      title = "MVNormal Dirichlet Process: R vs C++ Performance",
      subtitle = sprintf("Based on %d benchmark scenarios", nrow(bench_results) / 2)
    )

  return(list(
    dashboard = dashboard,
    speedup = p_speedup,
    scaling = p_scaling,
    dimension = p_dimension,
    clusters = p_clusters
  ))
}

#' Run complete MVNormal benchmark suite and generate report
#'
#' @export
run_mvnormal_benchmark_report <- function() {

  cat("\n")
  cat("========================================\n")
  cat("MVNormal Dirichlet Process Benchmarking\n")
  cat("========================================\n")

  # Check if C++ is available
  cpp_status <- get_cpp_status()
  if (!exists("conjugate_mvnormal_cluster_component_update_cpp")) {
    stop("MVNormal C++ implementation not available. Please recompile with C++ support.")
  }

  # 1. Quick benchmark
  cat("\n1. Running quick benchmark...\n")
  quick_results <- quick_mvnormal_benchmark()

  # 2. Component profiling
  cat("\n2. Profiling individual components...\n")
  component_results <- profile_mvnormal_components()

  # 3. Statistical equivalence
  cat("\n3. Testing statistical equivalence...\n")
  equiv_results <- test_mvnormal_statistical_equivalence()

  # 4. Comprehensive benchmarks (smaller scale for demo)
  cat("\n4. Running comprehensive benchmarks...\n")
  bench_results <- benchmark_mvnormal_comprehensive(
    n_obs_vec = c(100, 250, 500),
    dim_vec = c(2, 3, 5),
    n_iter_vec = c(100),
    n_clusters_vec = c(2, 3),
    n_reps = 3
  )

  # Generate summary statistics
  summary_stats <- bench_results %>%
    group_by(implementation) %>%
    summarise(
      mean_time = mean(time),
      median_time = median(time),
      mean_memory = mean(memory_mb, na.rm = TRUE),
      .groups = "drop"
    )

  # Calculate speedup by dimension
  speedup_by_dim <- bench_results %>%
    group_by(dimension, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop") %>%
    pivot_wider(names_from = implementation, values_from = mean_time) %>%
    mutate(speedup = R / `C++`)

  # Print summary report
  cat("\n")
  cat("========================================\n")
  cat("MVNORMAL BENCHMARK SUMMARY\n")
  cat("========================================\n")

  cat("\n1. OVERALL PERFORMANCE:\n")
  print(summary_stats)

  cat("\n2. SPEEDUP BY DIMENSION:\n")
  print(speedup_by_dim)

  cat("\n3. COMPONENT-LEVEL SPEEDUP:\n")
  comp_speedup <- component_results %>%
    pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
    mutate(speedup = R / `C++`)
  print(comp_speedup)

  cat("\n4. OVERALL SPEEDUP RANGE:\n")
  overall_speedup <- summary_stats$mean_time[summary_stats$implementation == "R"] /
    summary_stats$mean_time[summary_stats$implementation == "C++"]
  cat(sprintf("  Mean speedup: %.1fx\n", overall_speedup))
  cat(sprintf("  Min speedup: %.1fx\n", min(speedup_by_dim$speedup)))
  cat(sprintf("  Max speedup: %.1fx\n", max(speedup_by_dim$speedup)))

  # Create visualizations
  cat("\n5. Generating visualizations...\n")
  plots <- visualize_mvnormal_benchmarks(bench_results)

  # Display main dashboard
  print(plots$dashboard)

  # Save results
  results <- list(
    quick_benchmark = quick_results,
    component_profile = component_results,
    statistical_equivalence = equiv_results,
    comprehensive_benchmark = bench_results,
    summary_stats = summary_stats,
    speedup_by_dimension = speedup_by_dim,
    plots = plots
  )

  cat("\nBenchmarking complete!\n")

  return(invisible(results))
}

# Run the benchmark if sourced interactively
if (interactive()) {
  results <- run_mvnormal_benchmark_report()
}
