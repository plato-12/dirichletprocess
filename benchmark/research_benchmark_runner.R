# Research Benchmark Runner for Dirichlet Process Models
# ======================================================
#
# This script runs comprehensive research-quality benchmarks for all covariance models
# using the atime framework with publication-ready analysis and visualization.
#
# Features:
# - Comprehensive testing across multiple dimensions
# - Research-quality MCMC iterations and repetitions
# - Automatic result saving and report generation
# - Publication-ready plots using atime visualization
# - Detailed performance analysis and recommendations
#
# Usage:
#   source("benchmark/research_benchmark_runner.R")
#   run_full_research_benchmark()
#
# Author: Generated with Claude Code
# Date: 2025-01-16

# ==========================================
# SETUP AND DEPENDENCIES
# ==========================================

# Load required libraries
required_packages <- c("atime", "mvtnorm", "ggplot2", "dplyr", "data.table")

cat("=== RESEARCH BENCHMARK SETUP ===\n")
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
if (file.exists("R/benchmark_integration.R")) {
  devtools::load_all()
  source("benchmark/atime/benchmark-covariance-models-optimized.R")
} else {
  # If running from installed package
  library(dirichletprocess)
  if (file.exists("benchmark-covariance-models-optimized.R")) {
    source("benchmark-covariance-models-optimized.R")
  } else {
    stop("Cannot find benchmark system. Please run from package root directory.")
  }
}

cat("All packages loaded successfully.\n\n")

# ==========================================
# RESEARCH CONFIGURATION
# ==========================================

# Research-quality configuration
RESEARCH_CONFIG <- list(
  mcmc_iterations = 100,    # Thorough MCMC for research quality
  repetitions = 15,         # High statistical significance
  max_samples = 500         # Test scalability
)

# Publication configuration (even more thorough)
PUBLICATION_CONFIG <- list(
  mcmc_iterations = 200,    # Publication-quality MCMC
  repetitions = 20,         # Maximum statistical rigor
  max_samples = 1000        # Full scalability testing
)

# Quick research configuration (for testing)
QUICK_RESEARCH_CONFIG <- list(
  mcmc_iterations = 50,     # Reasonable for development
  repetitions = 10,         # Good statistical power
  max_samples = 200         # Moderate scalability
)

# ==========================================
# ENHANCED BENCHMARK FUNCTIONS
# ==========================================

#' Run comprehensive research benchmark with enhanced analysis
#'
#' @param config Configuration list (RESEARCH_CONFIG, PUBLICATION_CONFIG, etc.)
#' @param dimensions_list Vector of dimensions to test
#' @param save_results Whether to save results to files
#' @param generate_plots Whether to generate atime plots
#' @param output_dir Directory to save results (default: current directory)
#' @return List of benchmark results with analysis
#' @export
run_full_research_benchmark <- function(config = RESEARCH_CONFIG,
                                       dimensions_list = c(1, 2, 5, 10),
                                       save_results = TRUE,
                                       generate_plots = TRUE,
                                       output_dir = ".") {

  # Create output directory if needed
  if (save_results && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Print configuration
  cat("=== COMPREHENSIVE RESEARCH BENCHMARK ===\n")
  cat("Configuration:\n")
  cat("  MCMC iterations:", config$mcmc_iterations, "\n")
  cat("  Repetitions:", config$repetitions, "\n")
  cat("  Max samples:", config$max_samples, "\n")
  cat("  Dimensions to test:", paste(dimensions_list, collapse = ", "), "\n")
  cat("  Save results:", save_results, "\n")
  cat("  Generate plots:", generate_plots, "\n")
  cat("  Output directory:", output_dir, "\n\n")

  # Estimate time
  total_tests <- sum(ifelse(dimensions_list == 1, 3, 7))  # 3 for 1D, 7 for >1D
  estimated_minutes <- (total_tests * config$mcmc_iterations * config$repetitions) / 100
  cat("Estimated runtime:", round(estimated_minutes, 1), "minutes\n")
  cat("Starting benchmark...\n\n")

  # Start timing
  start_time <- Sys.time()

  # Run comprehensive benchmark
  benchmark_results <- run_comprehensive_benchmark(
    config = config,
    dimensions_list = dimensions_list
  )

  # Calculate total time
  end_time <- Sys.time()
  total_time <- end_time - start_time
  cat("Benchmark completed in:", format(total_time), "\n\n")

  # Generate enhanced analysis
  cat("=== GENERATING ENHANCED ANALYSIS ===\n")
  enhanced_results <- generate_enhanced_analysis(benchmark_results, config)

  # Generate and save plots if requested
  if (generate_plots) {
    cat("=== GENERATING ATIME PLOTS ===\n")
    plot_results <- generate_atime_plots(benchmark_results, output_dir, save_results)
    enhanced_results$plots <- plot_results
  }

  # Save results if requested
  if (save_results) {
    cat("=== SAVING RESULTS ===\n")
    save_research_results(enhanced_results, output_dir)
  }

  # Generate final report
  cat("=== GENERATING FINAL REPORT ===\n")
  generate_comprehensive_report(enhanced_results, total_time)

  return(enhanced_results)
}

#' Generate enhanced analysis with statistical insights
#'
#' @param benchmark_results Raw benchmark results
#' @param config Configuration used
#' @return Enhanced analysis with insights
generate_enhanced_analysis <- function(benchmark_results, config) {
  enhanced_analysis <- list(
    config = config,
    raw_results = benchmark_results,
    dimensional_analysis = list(),
    cross_dimensional_comparison = list(),
    statistical_summary = list(),
    performance_insights = list()
  )

  # Analyze each dimension
  for (dim_name in names(benchmark_results)) {
    cat("Analyzing", dim_name, "...\n")

    result <- benchmark_results[[dim_name]]
    analysis <- analyze_benchmark_results(result)

    # Add statistical analysis
    measurements <- result$measurements

    # Calculate confidence intervals
    confidence_intervals <- measurements[, .(
      mean_time = mean(median),
      sd_time = sd(median),
      ci_lower = mean(median) - 1.96 * sd(median) / sqrt(.N),
      ci_upper = mean(median) + 1.96 * sd(median) / sqrt(.N),
      n_measurements = .N
    ), by = expr.name]

    # Add to analysis
    analysis$confidence_intervals = confidence_intervals
    analysis$dimension = as.numeric(gsub("d", "", dim_name))

    enhanced_analysis$dimensional_analysis[[dim_name]] <- analysis
  }

  # Cross-dimensional comparison
  enhanced_analysis$cross_dimensional_comparison <- generate_cross_dimensional_analysis(
    enhanced_analysis$dimensional_analysis
  )

  # Generate insights
  enhanced_analysis$performance_insights <- generate_performance_insights(
    enhanced_analysis$dimensional_analysis
  )

  return(enhanced_analysis)
}

#' Generate cross-dimensional analysis
#'
#' @param dimensional_analysis Analysis results for each dimension
#' @return Cross-dimensional comparison insights
generate_cross_dimensional_analysis <- function(dimensional_analysis) {

  # Combine all results for cross-dimensional analysis
  all_results <- data.table()

  for (dim_name in names(dimensional_analysis)) {
    analysis <- dimensional_analysis[[dim_name]]
    dim_value <- analysis$dimension

    # Add dimension info to summary stats
    summary_with_dim <- analysis$summary_stats
    summary_with_dim$dimension <- dim_value
    summary_with_dim$dim_name <- dim_name

    all_results <- rbind(all_results, summary_with_dim, fill = TRUE)
  }

  # Analyze scaling behavior
  scaling_analysis <- list()

  # For each model, analyze how performance scales with dimension
  for (model in unique(all_results$expr.name)) {
    model_data <- all_results[expr.name == model]
    if (nrow(model_data) > 1) {
      # Simple linear model of time vs dimension
      lm_result <- lm(mean_time ~ dimension, data = model_data)
      scaling_analysis[[model]] <- list(
        slope = coef(lm_result)[2],
        intercept = coef(lm_result)[1],
        r_squared = summary(lm_result)$r.squared,
        scaling_interpretation = ifelse(coef(lm_result)[2] > 0.1, "Poor scaling",
                                       ifelse(coef(lm_result)[2] > 0.05, "Moderate scaling", "Good scaling"))
      )
    }
  }

  return(list(
    combined_results = all_results,
    scaling_analysis = scaling_analysis
  ))
}

#' Generate performance insights and recommendations
#'
#' @param dimensional_analysis Analysis results for each dimension
#' @return Performance insights and recommendations
generate_performance_insights <- function(dimensional_analysis) {
  insights <- list()

  # Find consistently fastest models across dimensions
  fastest_models <- sapply(dimensional_analysis, function(x) x$fastest_model)
  fastest_frequency <- table(fastest_models)

  # Find consistently slowest models
  slowest_models <- sapply(dimensional_analysis, function(x) x$slowest_model)
  slowest_frequency <- table(slowest_models)

  # Generate recommendations
  recommendations <- list(
    overall_fastest = names(fastest_frequency)[which.max(fastest_frequency)],
    overall_slowest = names(slowest_frequency)[which.max(slowest_frequency)],
    dimension_specific = list()
  )

  # Dimension-specific recommendations
  for (dim_name in names(dimensional_analysis)) {
    analysis <- dimensional_analysis[[dim_name]]
    dim_value <- analysis$dimension

    if (dim_value == 1) {
      recommendations$dimension_specific[[dim_name]] <- list(
        recommended = "V model for fastest performance, E model for balanced performance",
        avoid = "FULL model unless flexibility is critical",
        note = "Univariate models (E, V) significantly outperform FULL for 1D data"
      )
    } else {
      top_models <- head(analysis$summary_stats[order(mean_time)]$expr.name, 3)
      recommendations$dimension_specific[[dim_name]] <- list(
        recommended = paste("Top choices:", paste(top_models, collapse = ", ")),
        fastest = analysis$fastest_model,
        note = paste("For", dim_value, "dimensions, constrained models typically outperform FULL")
      )
    }
  }

  insights$recommendations <- recommendations
  insights$fastest_frequency <- fastest_frequency
  insights$slowest_frequency <- slowest_frequency

  return(insights)
}

#' Generate atime plots for visualization
#'
#' @param benchmark_results Raw benchmark results
#' @param output_dir Output directory for plots
#' @param save_plots Whether to save plots to files
#' @return List of plot objects
generate_atime_plots <- function(benchmark_results, output_dir = ".", save_plots = TRUE) {
  plots <- list()

  for (dim_name in names(benchmark_results)) {
    cat("Generating plots for", dim_name, "...\n")

    result <- benchmark_results[[dim_name]]

    # Generate atime plots - both time and memory
    tryCatch({
      # Time plot
      time_plot <- plot(result, unit = "seconds", log = "") +
        ggtitle(paste("Execution Time Comparison -", toupper(dim_name))) +
        theme_minimal() +
        theme(legend.position = "bottom")

      plots[[paste0(dim_name, "_time")]] <- time_plot

      # Memory plot
      memory_plot <- plot(result, unit = "kilobytes", log = "") +
        ggtitle(paste("Memory Usage Comparison -", toupper(dim_name))) +
        theme_minimal() +
        theme(legend.position = "bottom")

      plots[[paste0(dim_name, "_memory")]] <- memory_plot

      # Save plots if requested
      if (save_plots) {
        ggsave(
          filename = file.path(output_dir, paste0("time_comparison_", dim_name, ".png")),
          plot = time_plot,
          width = 10, height = 6, dpi = 300
        )

        ggsave(
          filename = file.path(output_dir, paste0("memory_comparison_", dim_name, ".png")),
          plot = memory_plot,
          width = 10, height = 6, dpi = 300
        )

        cat("  Saved plots for", dim_name, "\n")
      }

    }, error = function(e) {
      cat("  Warning: Could not generate plots for", dim_name, ":", e$message, "\n")
    })
  }

  # Generate combined scaling plot if multiple dimensions
  if (length(benchmark_results) > 1) {
    tryCatch({
      scaling_plot <- generate_scaling_plot(benchmark_results)
      plots$scaling <- scaling_plot

      if (save_plots) {
        ggsave(
          filename = file.path(output_dir, "performance_scaling.png"),
          plot = scaling_plot,
          width = 12, height = 8, dpi = 300
        )
        cat("  Saved scaling plot\n")
      }

    }, error = function(e) {
      cat("  Warning: Could not generate scaling plot:", e$message, "\n")
    })
  }

  return(plots)
}

#' Generate performance scaling plot across dimensions
#'
#' @param benchmark_results Raw benchmark results
#' @return ggplot object showing scaling behavior
generate_scaling_plot <- function(benchmark_results) {

  # Combine data from all dimensions
  combined_data <- data.table()

  for (dim_name in names(benchmark_results)) {
    measurements <- benchmark_results[[dim_name]]$measurements
    measurements$dimension <- as.numeric(gsub("d", "", dim_name))
    combined_data <- rbind(combined_data, measurements, fill = TRUE)
  }

  # Create scaling plot
  scaling_plot <- ggplot(combined_data, aes(x = dimension, y = median, color = expr.name)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
    scale_x_continuous(breaks = unique(combined_data$dimension)) +
    scale_y_log10() +
    labs(
      title = "Performance Scaling Across Dimensions",
      subtitle = "How execution time scales with data dimensionality",
      x = "Data Dimension",
      y = "Execution Time (seconds, log scale)",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )

  return(scaling_plot)
}

#' Save research results to files
#'
#' @param enhanced_results Enhanced analysis results
#' @param output_dir Output directory
save_research_results <- function(enhanced_results, output_dir) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Save R data file
  save(enhanced_results, file = file.path(output_dir, paste0("research_results_", timestamp, ".RData")))
  cat("Saved R data file\n")

  # Save CSV files for external analysis
  for (dim_name in names(enhanced_results$raw_results)) {
    measurements <- enhanced_results$raw_results[[dim_name]]$measurements
    write.csv(
      measurements,
      file = file.path(output_dir, paste0("measurements_", dim_name, "_", timestamp, ".csv")),
      row.names = FALSE
    )
  }
  cat("Saved CSV measurement files\n")

  # Save summary analysis
  summary_file <- file.path(output_dir, paste0("analysis_summary_", timestamp, ".txt"))
  sink(summary_file)
  generate_comprehensive_report(enhanced_results, NULL)
  sink()
  cat("Saved analysis summary\n")
}

#' Generate comprehensive research report
#'
#' @param enhanced_results Enhanced analysis results
#' @param total_time Total benchmark execution time
generate_comprehensive_report <- function(enhanced_results, total_time) {
  cat("=== COMPREHENSIVE RESEARCH BENCHMARK REPORT ===\n")
  cat("Generated:", format(Sys.time()), "\n")
  if (!is.null(total_time)) {
    cat("Total execution time:", format(total_time), "\n")
  }
  cat("Configuration:\n")
  cat("  MCMC iterations:", enhanced_results$config$mcmc_iterations, "\n")
  cat("  Repetitions:", enhanced_results$config$repetitions, "\n")
  cat("  Max samples:", enhanced_results$config$max_samples, "\n\n")

  # Dimensional analysis
  for (dim_name in names(enhanced_results$dimensional_analysis)) {
    analysis <- enhanced_results$dimensional_analysis[[dim_name]]

    cat(sprintf("=== %s ANALYSIS ===\n", toupper(dim_name)))
    cat("Models tested:", analysis$total_models_tested, "\n")
    cat("Sample sizes tested:", paste(analysis$sample_sizes_tested, collapse = ", "), "\n")
    cat("Fastest model:", analysis$fastest_model, "\n")
    cat("Slowest model:", analysis$slowest_model, "\n\n")

    cat("Performance Summary (with 95% confidence intervals):\n")
    ci_data <- analysis$confidence_intervals[order(mean_time)]
    for (i in 1:nrow(ci_data)) {
      row <- ci_data[i]
      cat(sprintf("  %s: %.3f ± %.3f seconds (95%% CI: %.3f - %.3f)\n",
                  row$expr.name, row$mean_time,
                  1.96 * row$sd_time / sqrt(row$n_measurements),
                  row$ci_lower, row$ci_upper))
    }
    cat("\n")
  }

  # Cross-dimensional insights
  if (length(enhanced_results$dimensional_analysis) > 1) {
    cat("=== CROSS-DIMENSIONAL INSIGHTS ===\n")
    scaling <- enhanced_results$cross_dimensional_comparison$scaling_analysis

    if (length(scaling) > 0) {
      cat("Scaling behavior (time vs dimension):\n")
      for (model in names(scaling)) {
        cat(sprintf("  %s: %s (R² = %.3f)\n",
                    model, scaling[[model]]$scaling_interpretation,
                    scaling[[model]]$r_squared))
      }
      cat("\n")
    }
  }

  # Recommendations
  cat("=== PERFORMANCE RECOMMENDATIONS ===\n")
  recommendations <- enhanced_results$performance_insights$recommendations

  cat("Overall fastest model:", recommendations$overall_fastest, "\n")
  cat("Overall slowest model:", recommendations$overall_slowest, "\n\n")

  cat("Dimension-specific recommendations:\n")
  for (dim_name in names(recommendations$dimension_specific)) {
    rec <- recommendations$dimension_specific[[dim_name]]
    cat(sprintf("  %s: %s\n", toupper(dim_name), rec$recommended))
    cat(sprintf("    Note: %s\n", rec$note))
  }

  cat("\n=== END REPORT ===\n")
}

# ==========================================
# QUICK EXECUTION FUNCTIONS
# ==========================================

#' Quick research benchmark (30-60 minutes)
#' @export
quick_research_benchmark <- function() {
  cat("Running quick research benchmark (estimated 30-60 minutes)...\n")
  return(run_full_research_benchmark(
    config = QUICK_RESEARCH_CONFIG,
    dimensions_list = c(1, 2, 5)
  ))
}

#' Standard research benchmark (1-3 hours)
#' @export
standard_research_benchmark <- function() {
  cat("Running standard research benchmark (estimated 1-3 hours)...\n")
  return(run_full_research_benchmark(
    config = RESEARCH_CONFIG,
    dimensions_list = c(1, 2, 5, 10)
  ))
}

#' Publication-quality benchmark (3-6 hours)
#' @export
publication_benchmark <- function() {
  cat("Running publication-quality benchmark (estimated 3-6 hours)...\n")
  return(run_full_research_benchmark(
    config = PUBLICATION_CONFIG,
    dimensions_list = c(1, 2, 5, 10, 20)
  ))
}

# ==========================================
# MAIN EXECUTION
# ==========================================

cat("=== RESEARCH BENCHMARK RUNNER LOADED ===\n")
cat("Available functions:\n")
cat("  quick_research_benchmark()    - Quick test (30-60 min)\n")
cat("  standard_research_benchmark() - Standard research (1-3 hours)\n")
cat("  publication_benchmark()       - Publication quality (3-6 hours)\n")
cat("  run_full_research_benchmark() - Custom configuration\n\n")

cat("Example usage:\n")
cat("  # Quick test\n")
cat("  results <- quick_research_benchmark()\n\n")
cat("  # Custom configuration\n")
cat("  results <- run_full_research_benchmark(\n")
cat("    config = RESEARCH_CONFIG,\n")
cat("    dimensions_list = c(1, 2, 5),\n")
cat("    save_results = TRUE,\n")
cat("    generate_plots = TRUE\n")
cat("  )\n\n")

cat("Ready to run benchmarks!\n")
