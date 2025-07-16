# Comprehensive MVNormal Covariance Models Benchmark Script
# =========================================================
#
# This script benchmarks all mvnormal covariance models for time and memory
# usage without using the atime package. It provides direct measurements
# and comprehensive analysis.
#
# Author: Generated with Claude Code
# Date: 2025-01-16

# ==========================================
# SETUP AND DEPENDENCIES
# ==========================================

# Load required libraries
required_packages <- c("mvtnorm", "ggplot2", "dplyr", "pryr", "microbenchmark")

cat("=== COMPREHENSIVE MVNORMAL BENCHMARK SETUP ===\n")
cat("Checking and installing required packages...\n")

for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Load dirichletprocess package
library(dirichletprocess)

# Try to load development version, fall back to installed version
if (file.exists("R/mvnormal_normal_wishart.R")) {
  devtools::load_all()
} else {
  library(dirichletprocess)
}

# Enable C++ for high-performance benchmarking
set_use_cpp(TRUE)

cat("All packages loaded successfully.\n\n")

# ==========================================
# BENCHMARK CONFIGURATION
# ==========================================

# Configuration for different benchmark intensities
BENCHMARK_CONFIGS <- list(
  QUICK = list(
    mcmc_iterations = 25,
    repetitions = 5,
    sample_sizes = c(50, 100),
    dimensions = c(1, 2, 5),
    warmup_runs = 2
  ),
  STANDARD = list(
    mcmc_iterations = 50,
    repetitions = 10,
    sample_sizes = c(50, 100, 200),
    dimensions = c(1, 2, 5, 10),
    warmup_runs = 3
  ),
  COMPREHENSIVE = list(
    mcmc_iterations = 100,
    repetitions = 15,
    sample_sizes = c(50, 100, 200, 500),
    dimensions = c(1, 2, 5, 10, 20),
    warmup_runs = 5
  )
)

# All covariance models to test
ALL_MODELS <- list(
  univariate = c("E", "V", "FULL"),
  multivariate = c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
)

# ==========================================
# UTILITY FUNCTIONS
# ==========================================

#' Generate benchmark data for testing
#'
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seed Random seed for reproducibility
#' @return Matrix of benchmark data
generate_benchmark_data <- function(n, d, seed = 42) {
  set.seed(seed)

  if (d == 1) {
    # Univariate case - mixture of normals
    data <- c(
      rnorm(n %/% 2, mean = -2, sd = 1),
      rnorm(n - n %/% 2, mean = 2, sd = 1)
    )
    return(matrix(data, ncol = 1))
  } else {
    # Multivariate case - create well-separated clusters
    k_clusters <- min(3, n %/% 10)
    cluster_sizes <- rep(n %/% k_clusters, k_clusters)
    cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)

    data_points <- c()
    for (i in 1:k_clusters) {
      mean_vec <- rep(0, d)
      mean_vec[1] <- (i - 2) * 3  # Separate along first dimension
      if (d > 1) mean_vec[2] <- (i - 2) * 2  # Separate along second dimension

      # Create covariance matrix
      sigma <- diag(d) * 0.8
      if (d > 1) {
        sigma[1, 2] <- 0.3  # Add some correlation
        sigma[2, 1] <- 0.3
      }

      # Generate cluster data
      cluster_data <- mvtnorm::rmvnorm(cluster_sizes[i], mean_vec, sigma)
      data_points <- rbind(data_points, cluster_data)
    }

    return(data_points[1:n, ])
  }
}

#' Create prior parameters for a given model and dimension
#'
#' @param dimensions Number of dimensions
#' @param model_name Covariance model name
#' @return List of prior parameters
create_prior_parameters <- function(dimensions, model_name) {
  # Validate model-dimension compatibility
  if (model_name %in% c("E", "V") && dimensions > 1) {
    stop(paste("Model", model_name, "is only for univariate data"))
  }

  if (model_name %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI") && dimensions == 1) {
    stop(paste("Model", model_name, "is only for multivariate data"))
  }

  # Create prior parameters
  prior_params <- list(
    mu0 = rep(0, dimensions),
    kappa0 = 1,
    nu = dimensions + 2,
    Lambda = diag(dimensions),
    covModel = model_name
  )

  return(prior_params)
}

#' Get valid models for a given dimension
#'
#' @param dimensions Number of dimensions
#' @return Vector of valid model names
get_valid_models <- function(dimensions) {
  if (dimensions == 1) {
    return(ALL_MODELS$univariate)
  } else {
    return(ALL_MODELS$multivariate)
  }
}


# ==========================================
# CORE BENCHMARK FUNCTIONS
# ==========================================

#' Run single benchmark iteration
#'
#' @param data_matrix Input data matrix
#' @param model_name Model name
#' @param mcmc_iterations Number of MCMC iterations
#' @param measure_memory Whether to measure memory usage
#' @return List with timing and memory results
run_single_benchmark <- function(data_matrix, model_name, mcmc_iterations, measure_memory = TRUE) {
  dimensions <- ncol(data_matrix)

  # Create prior parameters
  prior_params <- create_prior_parameters(dimensions, model_name)

  # Measure timing using system.time for more reliable results
  timing_result <- system.time({
    # Create mixing distribution and DP
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(data_matrix, md)
    dp <- Initialise(dp, numInitialClusters = 1)

    # Run MCMC
    result_obj <- Fit(dp, mcmc_iterations, progressBar = FALSE)
  })

  # Measure memory if requested
  memory_mb <- NA
  if (measure_memory && requireNamespace("pryr", quietly = TRUE)) {
    mem_before <- pryr::mem_used()

    # Run again for memory measurement
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(data_matrix, md)
    dp <- Initialise(dp, numInitialClusters = 1)
    result_obj <- Fit(dp, mcmc_iterations, progressBar = FALSE)

    mem_after <- pryr::mem_used()
    memory_mb <- as.numeric(mem_after - mem_before) / 1024^2
  }

  return(list(
    time_ms = (timing_result[["elapsed"]] * 1000),  # Convert to milliseconds
    memory_mb = memory_mb,
    n_clusters = result_obj$numberClusters,
    n_data_points = length(result_obj$clusterLabels),
    model = model_name,
    dimensions = dimensions,
    sample_size = nrow(data_matrix),
    mcmc_iterations = mcmc_iterations
  ))
}

#' Run benchmark for a specific configuration
#'
#' @param config Configuration list
#' @param model_name Model name
#' @param dimensions Number of dimensions
#' @param sample_size Sample size
#' @param verbose Whether to print progress
#' @return Data frame with benchmark results
run_model_benchmark <- function(config, model_name, dimensions, sample_size, verbose = TRUE) {

  if (verbose) {
    cat(sprintf("  Testing %s (d=%d, n=%d)... ", model_name, dimensions, sample_size))
  }

  # Generate data
  data_matrix <- generate_benchmark_data(sample_size, dimensions)

  # Run warmup iterations
  if (config$warmup_runs > 0) {
    for (i in 1:config$warmup_runs) {
      tryCatch({
        run_single_benchmark(data_matrix, model_name, 5, measure_memory = FALSE)
      }, error = function(e) {
        # Ignore warmup errors
      })
    }
  }

  # Run actual benchmark iterations
  results <- list()
  successful_runs <- 0

  for (rep in 1:config$repetitions) {
    tryCatch({
      result <- run_single_benchmark(data_matrix, model_name, config$mcmc_iterations, measure_memory = TRUE)
      result$repetition <- rep
      results[[length(results) + 1]] <- result
      successful_runs <- successful_runs + 1
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("Error in rep %d: %s ", rep, e$message))
      }
    })
  }

  if (verbose) {
    cat(sprintf("✓ %d/%d successful runs\n", successful_runs, config$repetitions))
  }

  if (length(results) == 0) {
    return(NULL)
  }

  # Convert to data frame
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(x, stringsAsFactors = FALSE)))

  return(results_df)
}

#' Run comprehensive benchmark across all configurations
#'
#' @param config Configuration list
#' @param verbose Whether to print progress
#' @return List with raw results and summary
run_comprehensive_benchmark <- function(config, verbose = TRUE) {

  if (verbose) {
    cat("=== COMPREHENSIVE BENCHMARK EXECUTION ===\n")
    cat("Configuration:\n")
    cat("  MCMC iterations:", config$mcmc_iterations, "\n")
    cat("  Repetitions:", config$repetitions, "\n")
    cat("  Sample sizes:", paste(config$sample_sizes, collapse = ", "), "\n")
    cat("  Dimensions:", paste(config$dimensions, collapse = ", "), "\n")
    cat("  Warmup runs:", config$warmup_runs, "\n\n")
  }

  # Calculate total combinations
  total_combinations <- 0
  for (d in config$dimensions) {
    valid_models <- get_valid_models(d)
    total_combinations <- total_combinations + length(valid_models) * length(config$sample_sizes)
  }

  estimated_minutes <- (total_combinations * config$repetitions * config$mcmc_iterations) / 1000

  if (verbose) {
    cat("Total combinations:", total_combinations, "\n")
    cat("Estimated time:", round(estimated_minutes, 1), "minutes\n\n")
  }

  # Run benchmarks
  all_results <- list()
  combination_count <- 0

  for (d in config$dimensions) {
    valid_models <- get_valid_models(d)

    if (verbose) {
      cat(sprintf("=== DIMENSION %d (%d models) ===\n", d, length(valid_models)))
    }

    for (model in valid_models) {
      for (sample_size in config$sample_sizes) {
        combination_count <- combination_count + 1

        result <- run_model_benchmark(config, model, d, sample_size, verbose)

        if (!is.null(result)) {
          all_results[[length(all_results) + 1]] <- result
        }

        if (verbose) {
          cat(sprintf("  Progress: %d/%d combinations completed\n", combination_count, total_combinations))
        }
      }
    }

    if (verbose) {
      cat("\n")
    }
  }

  # Combine all results
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
  } else {
    combined_results <- data.frame()
  }

  return(list(
    raw_results = combined_results,
    config = config,
    summary = generate_benchmark_summary(combined_results)
  ))
}

#' Generate summary statistics from benchmark results
#'
#' @param results_df Data frame with benchmark results
#' @return Summary statistics
generate_benchmark_summary <- function(results_df) {
  if (nrow(results_df) == 0) {
    return(list(message = "No successful benchmark results"))
  }

  # Summary by model and dimension
  summary_stats <- results_df %>%
    group_by(model, dimensions, sample_size) %>%
    summarise(
      n_runs = n(),
      mean_time_ms = mean(time_ms, na.rm = TRUE),
      median_time_ms = median(time_ms, na.rm = TRUE),
      sd_time_ms = sd(time_ms, na.rm = TRUE),
      min_time_ms = min(time_ms, na.rm = TRUE),
      max_time_ms = max(time_ms, na.rm = TRUE),
      mean_memory_mb = mean(memory_mb, na.rm = TRUE),
      median_memory_mb = median(memory_mb, na.rm = TRUE),
      mean_clusters = mean(n_clusters, na.rm = TRUE),
      .groups = "drop"
    )

  # Overall fastest and slowest
  fastest_overall <- summary_stats[which.min(summary_stats$median_time_ms), ]
  slowest_overall <- summary_stats[which.max(summary_stats$median_time_ms), ]

  # By dimension analysis
  by_dimension <- results_df %>%
    group_by(dimensions) %>%
    summarise(
      models_tested = n_distinct(model),
      fastest_model = model[which.min(time_ms)],
      slowest_model = model[which.max(time_ms)],
      fastest_time_ms = min(time_ms, na.rm = TRUE),
      slowest_time_ms = max(time_ms, na.rm = TRUE),
      .groups = "drop"
    )

  return(list(
    detailed_stats = summary_stats,
    fastest_overall = fastest_overall,
    slowest_overall = slowest_overall,
    by_dimension = by_dimension,
    total_successful_runs = nrow(results_df)
  ))
}

# ==========================================
# VISUALIZATION FUNCTIONS
# ==========================================

#' Create timing comparison plot
#'
#' @param results_df Data frame with benchmark results
#' @param title Plot title
#' @return ggplot object
create_timing_plot <- function(results_df, title = "Execution Time Comparison") {
  if (nrow(results_df) == 0) {
    return(ggplot() + ggtitle("No data available"))
  }

  # Calculate summary statistics for plotting
  plot_data <- results_df %>%
    group_by(model, dimensions, sample_size) %>%
    summarise(
      mean_time = mean(time_ms, na.rm = TRUE),
      median_time = median(time_ms, na.rm = TRUE),
      sd_time = sd(time_ms, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dimension_label = paste0("d=", dimensions),
      size_label = paste0("n=", sample_size)
    )

  ggplot(plot_data, aes(x = model, y = median_time, fill = model)) +
    geom_col(alpha = 0.7) +
    geom_errorbar(aes(ymin = median_time - sd_time, ymax = median_time + sd_time),
                  width = 0.2, alpha = 0.8) +
    facet_grid(dimension_label ~ size_label, scales = "free_y") +
    scale_y_log10() +
    labs(
      title = title,
      subtitle = "Median execution time with standard deviation error bars",
      x = "Covariance Model",
      y = "Execution Time (ms, log scale)",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
}

#' Create memory usage comparison plot
#'
#' @param results_df Data frame with benchmark results
#' @param title Plot title
#' @return ggplot object
create_memory_plot <- function(results_df, title = "Memory Usage Comparison") {
  if (nrow(results_df) == 0) {
    return(ggplot() + ggtitle("No data available"))
  }

  # Filter out NA memory values
  memory_data <- results_df %>%
    filter(!is.na(memory_mb)) %>%
    group_by(model, dimensions, sample_size) %>%
    summarise(
      mean_memory = mean(memory_mb, na.rm = TRUE),
      median_memory = median(memory_mb, na.rm = TRUE),
      sd_memory = sd(memory_mb, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dimension_label = paste0("d=", dimensions),
      size_label = paste0("n=", sample_size)
    )

  if (nrow(memory_data) == 0) {
    return(ggplot() + ggtitle("No memory data available"))
  }

  ggplot(memory_data, aes(x = model, y = median_memory, fill = model)) +
    geom_col(alpha = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, median_memory - sd_memory),
                      ymax = median_memory + sd_memory),
                  width = 0.2, alpha = 0.8) +
    facet_grid(dimension_label ~ size_label, scales = "free_y") +
    labs(
      title = title,
      subtitle = "Median memory usage with standard deviation error bars",
      x = "Covariance Model",
      y = "Memory Usage (MB)",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
}

#' Create scaling analysis plot
#'
#' @param results_df Data frame with benchmark results
#' @param title Plot title
#' @return ggplot object
create_scaling_plot <- function(results_df, title = "Performance Scaling Analysis") {
  if (nrow(results_df) == 0) {
    return(ggplot() + ggtitle("No data available"))
  }

  # Calculate median time for each configuration
  scaling_data <- results_df %>%
    group_by(model, dimensions, sample_size) %>%
    summarise(
      median_time = median(time_ms, na.rm = TRUE),
      .groups = "drop"
    )

  # Create separate plots for dimension and sample size scaling
  dim_scaling <- ggplot(scaling_data, aes(x = dimensions, y = median_time, color = model)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line(aes(group = model), alpha = 0.7) +
    facet_wrap(~paste0("n=", sample_size), scales = "free_y") +
    scale_y_log10() +
    labs(
      title = paste(title, "- Dimension Scaling"),
      x = "Dimensions",
      y = "Execution Time (ms, log scale)",
      color = "Model"
    ) +
    theme_minimal()

  return(dim_scaling)
}

#' Generate comprehensive visualization report
#'
#' @param benchmark_results Benchmark results list
#' @param save_plots Whether to save plots
#' @param output_dir Output directory
#' @return List of plot objects
generate_visualization_report <- function(benchmark_results, save_plots = TRUE, output_dir = ".") {
  results_df <- benchmark_results$raw_results

  plots <- list()

  # Timing plot
  plots$timing <- create_timing_plot(results_df, "MVNormal Covariance Models - Execution Time")

  # Memory plot
  plots$memory <- create_memory_plot(results_df, "MVNormal Covariance Models - Memory Usage")

  # Scaling plot
  plots$scaling <- create_scaling_plot(results_df, "MVNormal Covariance Models - Performance Scaling")

  # Save plots if requested
  if (save_plots) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    ggsave(
      filename = file.path(output_dir, paste0("timing_comparison_", timestamp, ".png")),
      plot = plots$timing,
      width = 12, height = 8, dpi = 300
    )

    ggsave(
      filename = file.path(output_dir, paste0("memory_comparison_", timestamp, ".png")),
      plot = plots$memory,
      width = 12, height = 8, dpi = 300
    )

    ggsave(
      filename = file.path(output_dir, paste0("scaling_analysis_", timestamp, ".png")),
      plot = plots$scaling,
      width = 12, height = 8, dpi = 300
    )

    cat("Plots saved to:", output_dir, "\n")
  }

  return(plots)
}

# ==========================================
# REPORT GENERATION
# ==========================================

#' Generate comprehensive text report
#'
#' @param benchmark_results Benchmark results list
#' @param execution_time Total execution time
generate_comprehensive_report <- function(benchmark_results, execution_time = NULL) {
  cat("=== COMPREHENSIVE MVNORMAL BENCHMARK REPORT ===\n")
  cat("Generated:", format(Sys.time()), "\n")

  if (!is.null(execution_time)) {
    cat("Total execution time:", format(execution_time), "\n")
  }

  config <- benchmark_results$config
  cat("Configuration:\n")
  cat("  MCMC iterations:", config$mcmc_iterations, "\n")
  cat("  Repetitions:", config$repetitions, "\n")
  cat("  Sample sizes:", paste(config$sample_sizes, collapse = ", "), "\n")
  cat("  Dimensions:", paste(config$dimensions, collapse = ", "), "\n")
  cat("  Warmup runs:", config$warmup_runs, "\n\n")

  summary <- benchmark_results$summary

  if (is.null(summary$detailed_stats) || nrow(summary$detailed_stats) == 0) {
    cat("No successful benchmark results to report.\n")
    return()
  }

  cat("OVERALL RESULTS:\n")
  cat("  Total successful runs:", summary$total_successful_runs, "\n")
  cat("  Fastest model overall:", summary$fastest_overall$model,
      sprintf("(%.2f ms)", summary$fastest_overall$median_time_ms), "\n")
  cat("  Slowest model overall:", summary$slowest_overall$model,
      sprintf("(%.2f ms)", summary$slowest_overall$median_time_ms), "\n\n")

  # Results by dimension
  cat("RESULTS BY DIMENSION:\n")
  for (i in 1:nrow(summary$by_dimension)) {
    dim_info <- summary$by_dimension[i, ]
    cat(sprintf("  %dD: %d models tested\n", dim_info$dimensions, dim_info$models_tested))
    cat(sprintf("    Fastest: %s (%.2f ms)\n", dim_info$fastest_model, dim_info$fastest_time_ms))
    cat(sprintf("    Slowest: %s (%.2f ms)\n", dim_info$slowest_model, dim_info$slowest_time_ms))
  }
  cat("\n")

  # Detailed performance statistics
  cat("DETAILED PERFORMANCE STATISTICS:\n")
  detailed_stats <- summary$detailed_stats[order(summary$detailed_stats$median_time_ms), ]

  for (i in 1:nrow(detailed_stats)) {
    stat <- detailed_stats[i, ]
    cat(sprintf("  %s (d=%d, n=%d): %.2f ± %.2f ms (%.2f MB memory)\n",
                stat$model, stat$dimensions, stat$sample_size,
                stat$median_time_ms, stat$sd_time_ms, stat$mean_memory_mb))
  }

  cat("\n=== END REPORT ===\n")
}

# ==========================================
# MAIN EXECUTION FUNCTIONS
# ==========================================

#' Run quick benchmark (5-10 minutes)
#'
#' @param save_results Whether to save results
#' @param generate_plots Whether to generate plots
#' @return Benchmark results
run_quick_benchmark <- function(save_results = TRUE, generate_plots = TRUE) {
  cat("Running quick benchmark (estimated 5-10 minutes)...\n")

  start_time <- Sys.time()
  results <- run_comprehensive_benchmark(BENCHMARK_CONFIGS$QUICK, verbose = TRUE)
  end_time <- Sys.time()

  execution_time <- end_time - start_time

  if (generate_plots) {
    plots <- generate_visualization_report(results, save_plots = save_results)
    results$plots <- plots
  }

  generate_comprehensive_report(results, execution_time)

  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0("benchmark_results_", timestamp, ".RData"))
    cat("Results saved to: benchmark_results_", timestamp, ".RData\n")
  }

  return(results)
}

#' Run standard benchmark (15-30 minutes)
#'
#' @param save_results Whether to save results
#' @param generate_plots Whether to generate plots
#' @return Benchmark results
run_standard_benchmark <- function(save_results = TRUE, generate_plots = TRUE) {
  cat("Running standard benchmark (estimated 15-30 minutes)...\n")

  start_time <- Sys.time()
  results <- run_comprehensive_benchmark(BENCHMARK_CONFIGS$STANDARD, verbose = TRUE)
  end_time <- Sys.time()

  execution_time <- end_time - start_time

  if (generate_plots) {
    plots <- generate_visualization_report(results, save_plots = save_results)
    results$plots <- plots
  }

  generate_comprehensive_report(results, execution_time)

  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0("benchmark_results_", timestamp, ".RData"))
    cat("Results saved to: benchmark_results_", timestamp, ".RData\n")
  }

  return(results)
}

#' Run comprehensive benchmark (30-60 minutes)
#'
#' @param save_results Whether to save results
#' @param generate_plots Whether to generate plots
#' @return Benchmark results
run_comprehensive_benchmark_full <- function(save_results = TRUE, generate_plots = TRUE) {
  cat("Running comprehensive benchmark (estimated 30-60 minutes)...\n")

  start_time <- Sys.time()
  results <- run_comprehensive_benchmark(BENCHMARK_CONFIGS$COMPREHENSIVE, verbose = TRUE)
  end_time <- Sys.time()

  execution_time <- end_time - start_time

  if (generate_plots) {
    plots <- generate_visualization_report(results, save_plots = save_results)
    results$plots <- plots
  }

  generate_comprehensive_report(results, execution_time)

  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0("benchmark_results_", timestamp, ".RData"))
    cat("Results saved to: benchmark_results_", timestamp, ".RData\n")
  }

  return(results)
}

# ==========================================
# SCRIPT INITIALIZATION
# ==========================================

cat("=== COMPREHENSIVE MVNORMAL BENCHMARK LOADED ===\n")
cat("Available functions:\n")
cat("  run_quick_benchmark()        - Quick test (5-10 min)\n")
cat("  run_standard_benchmark()     - Standard test (15-30 min)\n")
cat("  run_comprehensive_benchmark_full() - Full test (30-60 min)\n")
cat("  run_comprehensive_benchmark() - Custom configuration\n\n")

cat("Example usage:\n")
cat("  # Quick benchmark\n")
cat("  results <- run_quick_benchmark()\n\n")
cat("  # Standard benchmark\n")
cat("  results <- run_standard_benchmark()\n\n")
cat("  # Custom configuration\n")
cat("  results <- run_comprehensive_benchmark(BENCHMARK_CONFIGS$QUICK)\n\n")

cat("Models that will be tested:\n")
cat("  1D (univariate):", paste(ALL_MODELS$univariate, collapse = ", "), "\n")
cat("  2D+ (multivariate):", paste(ALL_MODELS$multivariate, collapse = ", "), "\n\n")

cat("Ready to run benchmarks!\n")
