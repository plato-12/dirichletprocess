# benchmark/visualize_performance.R

create_performance_report <- function(results) {
  # Convert results to data frame
  perf_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      distribution = x$distribution,
      sample_size = x$sample_size,
      iterations = x$iterations,
      speedup = x$speedup
    )
  }))

  # Speedup by sample size
  p1 <- ggplot(perf_df, aes(x = sample_size, y = speedup, color = distribution)) +
    geom_line() +
    geom_point() +
    facet_wrap(~iterations, scales = "free_y") +
    labs(title = "C++ Speedup by Sample Size and Iterations",
         x = "Sample Size", y = "Speedup Factor") +
    theme_minimal()

  ggsave("performance_speedup_by_size.png", p1, width = 12, height = 8)

  # Average speedup by distribution
  avg_speedup <- aggregate(speedup ~ distribution, perf_df, mean)

  p2 <- ggplot(avg_speedup, aes(x = reorder(distribution, speedup), y = speedup)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Average C++ Speedup by Distribution",
         x = "Distribution", y = "Average Speedup Factor") +
    theme_minimal()

  ggsave("performance_speedup_by_distribution.png", p2, width = 8, height = 6)

  # Create summary report
  cat("\n=== PERFORMANCE SUMMARY ===\n")
  print(avg_speedup)
  cat("\nOverall average speedup:", round(mean(perf_df$speedup), 2), "x\n")
}

plot_scaling_results <- function(scaling_results) {
  # Combine results
  scaling_df <- do.call(rbind, lapply(names(scaling_results), function(dist) {
    df <- scaling_results[[dist]]
    df$distribution <- dist
    df
  }))

  # Time scaling plot
  p1 <- ggplot(scaling_df, aes(x = n)) +
    geom_line(aes(y = r_time, color = "R"), size = 1.2) +
    geom_line(aes(y = cpp_time, color = "C++"), size = 1.2) +
    facet_wrap(~distribution) +
    scale_y_log10() +
    labs(title = "Computation Time Scaling",
         x = "Sample Size", y = "Time (seconds, log scale)") +
    theme_minimal()

  ggsave("performance_scaling.png", p1, width = 10, height = 6)

  # Speedup scaling plot
  p2 <- ggplot(scaling_df, aes(x = n, y = speedup, color = distribution)) +
    geom_line(size = 1.2) +
    geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
    labs(title = "Speedup Factor vs Sample Size",
         x = "Sample Size", y = "Speedup Factor") +
    theme_minimal()

  ggsave("performance_speedup_scaling.png", p2, width = 8, height = 6)
}

# Create comprehensive performance dashboard
create_performance_dashboard <- function(perf_results, scaling_results, memory_results) {
  # Create output directory
  dir.create("performance_dashboard", showWarnings = FALSE)

  # 1. Distribution comparison heatmap
  perf_df <- do.call(rbind, lapply(perf_results, function(x) {
    data.frame(
      distribution = x$distribution,
      sample_size = x$sample_size,
      iterations = x$iterations,
      speedup = x$speedup
    )
  }))

  p_heatmap <- ggplot(perf_df, aes(x = factor(sample_size),
                                   y = factor(iterations),
                                   fill = speedup)) +
    geom_tile() +
    facet_wrap(~distribution) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 1, name = "Speedup") +
    labs(title = "Performance Heatmap: C++ Speedup Factors",
         x = "Sample Size", y = "Iterations") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave("performance_dashboard/speedup_heatmap.png", p_heatmap, width = 12, height = 8)

  # 2. Box plot of speedup distributions
  p_boxplot <- ggplot(perf_df, aes(x = distribution, y = speedup)) +
    geom_boxplot(fill = "lightblue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Distribution of C++ Speedup Factors",
         x = "Distribution", y = "Speedup Factor") +
    theme_minimal()

  ggsave("performance_dashboard/speedup_boxplot.png", p_boxplot, width = 10, height = 6)

  # 3. Performance by sample size (log scale)
  p_logscale <- ggplot(perf_df, aes(x = sample_size, y = speedup,
                                    color = distribution)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    scale_x_log10() +
    labs(title = "Speedup vs Sample Size (Log Scale)",
         x = "Sample Size (log scale)", y = "Speedup Factor") +
    theme_minimal()

  ggsave("performance_dashboard/speedup_logscale.png", p_logscale, width = 10, height = 6)

  # 4. Memory efficiency plot (if available)
  if (!is.null(memory_results)) {
    mem_df <- data.frame(
      distribution = names(memory_results),
      memory_ratio = sapply(memory_results, `[[`, "memory_ratio")
    )

    p_memory <- ggplot(mem_df, aes(x = distribution, y = memory_ratio)) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(title = "Memory Efficiency: R/C++ Memory Usage Ratio",
           x = "Distribution", y = "Memory Ratio (R/C++)") +
      theme_minimal()

    ggsave("performance_dashboard/memory_efficiency.png", p_memory, width = 8, height = 6)
  }

  # 5. Generate HTML report
  generate_html_performance_report(perf_df, scaling_results, memory_results)
}

# Generate HTML performance report
generate_html_performance_report <- function(perf_df, scaling_results, memory_results) {
  # Build HTML content using paste0
  html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
    <title>Dirichlet Process C++ Performance Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #333; }
        h2 { color: #666; }
        .metric {
            background-color: #f0f0f0;
            padding: 20px;
            margin: 10px 0;
            border-radius: 5px;
        }
        .metric-value {
            font-size: 2em;
            font-weight: bold;
            color: #0066cc;
        }
        img { max-width: 100%; height: auto; margin: 20px 0; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>Dirichlet Process C++ Performance Report</h1>
    <p>Generated: ', format(Sys.Date(), "%Y-%m-%d"), '</p>

    <h2>Executive Summary</h2>
    <div class="metric">
        <p>Overall Average Speedup</p>
        <p class="metric-value">', sprintf("%.2fx", mean(perf_df$speedup)), '</p>
    </div>
')

  # Add distribution-specific metrics
  for (dist in unique(perf_df$distribution)) {
    dist_speedup <- mean(perf_df$speedup[perf_df$distribution == dist])
    html_content <- paste0(html_content, '
    <div class="metric">
        <p>', dist, ' Distribution Average Speedup</p>
        <p class="metric-value">', sprintf("%.2fx", dist_speedup), '</p>
    </div>
')
  }

  html_content <- paste0(html_content, '
    <h2>Performance Visualizations</h2>
    <img src="speedup_heatmap.png" alt="Speedup Heatmap">
    <img src="speedup_boxplot.png" alt="Speedup Distribution">
    <img src="speedup_logscale.png" alt="Speedup vs Sample Size">
')

  if (!is.null(memory_results)) {
    html_content <- paste0(html_content, '
    <img src="memory_efficiency.png" alt="Memory Efficiency">
')
  }

  html_content <- paste0(html_content, '
    <h2>Detailed Results Table</h2>
    <table>
        <tr>
            <th>Distribution</th>
            <th>Sample Size</th>
            <th>Iterations</th>
            <th>Speedup Factor</th>
        </tr>
')

  # Add table rows
  for (i in 1:nrow(perf_df)) {
    html_content <- paste0(html_content, sprintf('
        <tr>
            <td>%s</td>
            <td>%d</td>
            <td>%d</td>
            <td>%.2fx</td>
        </tr>',
                                              perf_df$distribution[i],
                                              perf_df$sample_size[i],
                                              perf_df$iterations[i],
                                              perf_df$speedup[i]
    ))
  }

  html_content <- paste0(html_content, '
    </table>
</body>
</html>
')

  writeLines(html_content, "performance_dashboard/performance_report.html")
  cat("HTML report generated: performance_dashboard/performance_report.html\n")
}

# Helper to create timing comparison plots
plot_timing_comparison <- function(benchmark_results) {
  # Extract timing data
  timing_df <- do.call(rbind, lapply(names(benchmark_results), function(dist) {
    bench <- benchmark_results[[dist]]$benchmark
    data.frame(
      distribution = dist,
      implementation = bench$expr,
      time = bench$time / 1e9  # Convert to seconds
    )
  }))

  # Create violin plot
  p <- ggplot(timing_df, aes(x = distribution, y = time, fill = implementation)) +
    geom_violin(alpha = 0.7) +
    scale_y_log10() +
    labs(title = "Timing Distribution Comparison",
         x = "Distribution", y = "Time (seconds, log scale)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
