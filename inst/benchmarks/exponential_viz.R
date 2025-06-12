# Visualization functions for Exponential DP benchmarks

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

#' Create comprehensive visualization dashboard
#'
#' @param bench_data Benchmark results from benchmark_exponential_comprehensive
#' @param memory_data Memory profiling results
#' @param component_data Component benchmark results
#' @return List of ggplot objects
create_benchmark_dashboard <- function(bench_data, memory_data, component_data) {

  # Set consistent theme
  theme_set(theme_minimal(base_size = 12) +
              theme(
                legend.position = "bottom",
                panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 14)
              ))

  # Color palette
  impl_colors <- c("R" = "#E74C3C", "C++" = "#3498DB")

  # 1. Time vs Dataset Size (log-log plot)
  p_time_scaling <- bench_data %>%
    group_by(n_obs, n_iter, implementation) %>%
    summarise(
      mean_time = mean(time),
      se_time = sd(time) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = n_obs, y = mean_time, color = implementation)) +
    geom_line(aes(linetype = as.factor(n_iter)), size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                  width = 0.05, alpha = 0.5) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = function(x) sprintf("%.2f", x)) +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Execution Time Scaling",
      x = "Number of Observations",
      y = "Time (seconds)",
      color = "Implementation",
      linetype = "Iterations"
    ) +
    annotation_logticks()

  # 2. Speedup heatmap
  speedup_data <- bench_data %>%
    group_by(n_obs, n_iter, n_true_clusters, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop") %>%
    pivot_wider(names_from = implementation, values_from = mean_time) %>%
    mutate(speedup = R / `C++`)

  p_speedup_heatmap <- ggplot(speedup_data,
                              aes(x = as.factor(n_obs), y = as.factor(n_iter),
                                  fill = speedup)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1fx", speedup)), size = 3) +
    facet_wrap(~ n_true_clusters, labeller = label_both) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 1, limits = c(0.5, NA)) +
    labs(
      title = "C++ Speedup Factor",
      x = "Number of Observations",
      y = "Iterations",
      fill = "Speedup"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # 3. Memory usage comparison
  p_memory <- bench_data %>%
    group_by(n_obs, implementation) %>%
    summarise(
      mean_memory = mean(memory_mb),
      se_memory = sd(memory_mb) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = n_obs, y = mean_memory, fill = implementation)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_memory - se_memory,
                      ymax = mean_memory + se_memory),
                  position = position_dodge(0.8), width = 0.3) +
    scale_x_continuous(breaks = unique(bench_data$n_obs)) +
    scale_fill_manual(values = impl_colors) +
    labs(
      title = "Memory Usage",
      x = "Number of Observations",
      y = "Memory (MB)",
      fill = "Implementation"
    )

  # 4. Component-level performance
  p_components <- component_data %>%
    mutate(component = factor(component,
                              levels = c("Likelihood", "ClusterAssignment", "ParameterUpdate"))) %>%
    ggplot(aes(x = component, y = mean_time_ms, fill = implementation)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_time_ms - sd_time_ms,
                      ymax = mean_time_ms + sd_time_ms),
                  position = position_dodge(0.8), width = 0.3) +
    scale_fill_manual(values = impl_colors) +
    labs(
      title = "Component-Level Performance",
      x = "Component",
      y = "Time (milliseconds)",
      fill = "Implementation"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # 5. Clustering accuracy
  p_clustering <- bench_data %>%
    group_by(n_true_clusters, n_obs, implementation) %>%
    summarise(
      mean_found = mean(clusters_found),
      se_found = sd(clusters_found) / sqrt(n()),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = n_true_clusters, y = mean_found, color = implementation)) +
    geom_point(aes(size = n_obs), position = position_dodge(0.3), alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_found - se_found, ymax = mean_found + se_found),
                  position = position_dodge(0.3), width = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = impl_colors) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = "Clustering Accuracy",
      x = "True Number of Clusters",
      y = "Mean Clusters Found",
      color = "Implementation",
      size = "N Obs"
    )

  # 6. Memory timeline (if available)
  if (!is.null(memory_data$r_memory$during_fit)) {
    mem_timeline <- data.frame(
      iteration = rep(seq(10, 100, 10), 2),
      memory_mb = c(memory_data$r_memory$during_fit / 1024^2,
                    memory_data$cpp_memory$during_fit / 1024^2),
      implementation = rep(c("R", "C++"), each = 10)
    )

    p_mem_timeline <- ggplot(mem_timeline,
                             aes(x = iteration, y = memory_mb, color = implementation)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      scale_color_manual(values = impl_colors) +
      labs(
        title = "Memory Usage During Fitting",
        x = "Iteration Progress (%)",
        y = "Memory (MB)",
        color = "Implementation"
      )
  } else {
    p_mem_timeline <- NULL
  }

  # Create combined dashboard
  dashboard <- (p_time_scaling + p_speedup_heatmap) /
    (p_memory + p_components) /
    (p_clustering + (p_mem_timeline %||% plot_spacer()))

  dashboard <- dashboard +
    plot_annotation(
      title = "Exponential Distribution Performance Benchmarks",
      subtitle = "Comparing R vs C++ Implementations",
      theme = theme(plot.title = element_text(size = 18, face = "bold"))
    )

  return(list(
    dashboard = dashboard,
    time_scaling = p_time_scaling,
    speedup_heatmap = p_speedup_heatmap,
    memory = p_memory,
    components = p_components,
    clustering = p_clustering,
    memory_timeline = p_mem_timeline
  ))
}

#' Create detailed performance report plot
#'
#' @param bench_data Benchmark results
#' @return ggplot object
create_performance_report <- function(bench_data) {

  # Calculate summary statistics
  summary_stats <- bench_data %>%
    group_by(implementation) %>%
    summarise(
      mean_time = mean(time),
      median_time = median(time),
      q25_time = quantile(time, 0.25),
      q75_time = quantile(time, 0.75),
      mean_memory = mean(memory_mb),
      median_memory = median(memory_mb),
      .groups = "drop"
    )

  # Create text summary
  speedup_stats <- bench_data %>%
    select(scenario_id, rep, implementation, time) %>%
    pivot_wider(names_from = implementation, values_from = time) %>%
    summarise(
      mean_speedup = mean(R / `C++`, na.rm = TRUE),
      median_speedup = median(R / `C++`, na.rm = TRUE),
      min_speedup = min(R / `C++`, na.rm = TRUE),
      max_speedup = max(R / `C++`, na.rm = TRUE)
    )

  # Create a text plot with summary
  p_summary <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.9,
             label = "Performance Summary",
             size = 6, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.75,
             label = sprintf("Mean Speedup: %.1fx", speedup_stats$mean_speedup),
             size = 5) +
    annotate("text", x = 0.5, y = 0.65,
             label = sprintf("Median Speedup: %.1fx", speedup_stats$median_speedup),
             size = 5) +
    annotate("text", x = 0.5, y = 0.55,
             label = sprintf("Range: %.1fx - %.1fx",
                             speedup_stats$min_speedup,
                             speedup_stats$max_speedup),
             size = 5) +
    annotate("text", x = 0.5, y = 0.4,
             label = "Memory Reduction",
             size = 5, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.3,
             label = sprintf("R: %.1f MB (median)",
                             summary_stats$median_memory[summary_stats$implementation == "R"]),
             size = 4) +
    annotate("text", x = 0.5, y = 0.2,
             label = sprintf("C++: %.1f MB (median)",
                             summary_stats$median_memory[summary_stats$implementation == "C++"]),
             size = 4) +
    xlim(0, 1) + ylim(0, 1)

  return(p_summary)
}

#' Save all plots to files
#'
#' @param plots List of plots from create_benchmark_dashboard
#' @param output_dir Directory to save plots
save_benchmark_plots <- function(plots, output_dir = "benchmark_results/exponential") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save dashboard
  ggsave(
    file.path(output_dir, "dashboard.png"),
    plots$dashboard,
    width = 16, height = 20, dpi = 300
  )

  # Save individual plots
  plot_names <- c("time_scaling", "speedup_heatmap", "memory",
                  "components", "clustering", "memory_timeline")

  for (name in plot_names) {
    if (!is.null(plots[[name]])) {
      ggsave(
        file.path(output_dir, paste0(name, ".png")),
        plots[[name]],
        width = 8, height = 6, dpi = 300
      )
    }
  }

  cat("Plots saved to:", output_dir, "\n")
}

# Example usage function
run_exponential_visualization <- function() {
  # First run the benchmarks (from previous script)
  source("benchmark_exponential.R")
  results <- run_and_report()

  # Create visualizations
  plots <- create_benchmark_dashboard(
    results$benchmarks,
    results$memory,
    results$components
  )

  # Save plots
  save_benchmark_plots(plots)

  # Display dashboard
  print(plots$dashboard)

  return(plots)
}
