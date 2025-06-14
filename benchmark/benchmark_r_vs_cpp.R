# Comprehensive benchmark comparing R and C++ backends

library(dirichletprocess)
library(microbenchmark)
library(ggplot2)
library(profmem)
library(dplyr)

# Benchmark Configuration
BENCHMARK_CONFIG <- list(
  data_sizes = c(50, 100, 200, 500, 1000),
  n_iterations = c(100, 200, 500),
  n_replicates = 5,
  random_seed = 42
)

# Helper Functions
create_benchmark_data <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Create data with clear clustering structure
  cluster1 <- rnorm(n %/% 3, mean = -2, sd = 0.5)
  cluster2 <- rnorm(n %/% 3, mean = 0, sd = 0.7)
  cluster3 <- rnorm(n - 2 * (n %/% 3), mean = 2, sd = 0.6)
  c(cluster1, cluster2, cluster3)
}

setup_dirichlet_process <- function(data, backend = "r") {
  set_use_cpp(backend == "cpp")
  DirichletProcessGaussian(data)
}

# Memory Usage Benchmark
benchmark_memory <- function(data_size, n_iter, backend) {
  data <- create_benchmark_data(data_size, BENCHMARK_CONFIG$random_seed)

  if (backend == "r") {
    set_use_cpp(FALSE)
    prof <- profmem({
      dp <- DirichletProcessGaussian(data)
      dp <- Fit(dp, n_iter, progressBar = FALSE)
    })
  } else {
    set_use_cpp(TRUE)
    prof <- profmem({
      dp <- DirichletProcessGaussian(data)
      dp <- Fit(dp, n_iter, progressBar = FALSE)
    })
  }

  total_memory <- sum(prof$bytes, na.rm = TRUE)
  return(list(
    backend = backend,
    data_size = data_size,
    n_iter = n_iter,
    memory_bytes = total_memory,
    memory_mb = total_memory / (1024^2)
  ))
}

# Speed Benchmark
benchmark_speed <- function(data_size, n_iter, times = 3) {
  data <- create_benchmark_data(data_size, BENCHMARK_CONFIG$random_seed)

  # Benchmark R backend
  mb_r <- microbenchmark(
    r_backend = {
      set_use_cpp(FALSE)
      dp <- DirichletProcessGaussian(data)
      dp <- Fit(dp, n_iter, progressBar = FALSE)
    },
    times = times
  )

  # Benchmark C++ backend
  mb_cpp <- microbenchmark(
    cpp_backend = {
      set_use_cpp(TRUE)
      dp <- DirichletProcessGaussian(data)
      dp <- Fit(dp, n_iter, progressBar = FALSE)
    },
    times = times
  )

  return(list(
    data_size = data_size,
    n_iter = n_iter,
    r_median_ms = median(mb_r$time) / 1e6,
    cpp_median_ms = median(mb_cpp$time) / 1e6,
    r_times = mb_r$time / 1e6,
    cpp_times = mb_cpp$time / 1e6,
    speedup = median(mb_r$time) / median(mb_cpp$time)
  ))
}

# Quality Comparison
compare_quality <- function(data_size, n_iter, seed = 123) {
  data <- create_benchmark_data(data_size, seed)

  # R backend
  set.seed(seed)
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(data)
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  # C++ backend
  set.seed(seed)
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(data)
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  return(list(
    data_size = data_size,
    n_iter = n_iter,
    r_clusters = dp_r$numberClusters,
    cpp_clusters = dp_cpp$numberClusters,
    r_alpha = dp_r$alpha,
    cpp_alpha = dp_cpp$alpha,
    r_loglik = sum(log(PredictiveArray(dp_r, dp_r$data))),
    cpp_loglik = sum(log(PredictiveArray(dp_cpp, dp_cpp$data)))
  ))
}

# Main Benchmark Function
run_comprehensive_benchmark <- function() {
  cat("Starting Comprehensive R vs C++ Benchmark\n")
  cat("==========================================\n\n")

  # Initialize results
  speed_results <- list()
  memory_results <- list()
  quality_results <- list()

  total_combinations <- length(BENCHMARK_CONFIG$data_sizes) * length(BENCHMARK_CONFIG$n_iterations)
  current_combination <- 0

  for (data_size in BENCHMARK_CONFIG$data_sizes) {
    for (n_iter in BENCHMARK_CONFIG$n_iterations) {
      current_combination <- current_combination + 1
      cat(sprintf("Progress: %d/%d - Data size: %d, Iterations: %d\n",
                  current_combination, total_combinations, data_size, n_iter))

      # Speed benchmark
      tryCatch({
        speed_result <- benchmark_speed(data_size, n_iter, BENCHMARK_CONFIG$n_replicates)
        speed_results[[length(speed_results) + 1]] <- speed_result
      }, error = function(e) {
        cat("Speed benchmark failed:", e$message, "\n")
      })

      # Memory benchmark
      tryCatch({
        memory_r <- benchmark_memory(data_size, n_iter, "r")
        memory_cpp <- benchmark_memory(data_size, n_iter, "cpp")
        memory_results[[length(memory_results) + 1]] <- memory_r
        memory_results[[length(memory_results) + 1]] <- memory_cpp
      }, error = function(e) {
        cat("Memory benchmark failed:", e$message, "\n")
      })

      # Quality comparison
      tryCatch({
        quality_result <- compare_quality(data_size, n_iter)
        quality_results[[length(quality_results) + 1]] <- quality_result
      }, error = function(e) {
        cat("Quality comparison failed:", e$message, "\n")
      })
    }
  }

  return(list(
    speed = speed_results,
    memory = memory_results,
    quality = quality_results
  ))
}

# Results Analysis and Visualization
analyze_benchmark_results <- function(results) {
  cat("\nBenchmark Results Analysis\n")
  cat("=========================\n\n")

  # Process speed results
  if (length(results$speed) > 0) {
    speed_df <- do.call(rbind, lapply(results$speed, function(x) {
      data.frame(
        data_size = x$data_size,
        n_iter = x$n_iter,
        r_time = x$r_median_ms,
        cpp_time = x$cpp_median_ms,
        speedup = x$speedup
      )
    }))

    cat("Speed Results Summary:\n")
    cat("----------------------\n")
    print(speed_df)
    cat(sprintf("Average speedup: %.2fx\n", mean(speed_df$speedup)))
    cat(sprintf("Maximum speedup: %.2fx\n", max(speed_df$speedup)))
    cat(sprintf("Minimum speedup: %.2fx\n", min(speed_df$speedup)))
  }

  # Process memory results
  if (length(results$memory) > 0) {
    memory_df <- do.call(rbind, lapply(results$memory, function(x) {
      data.frame(
        data_size = x$data_size,
        n_iter = x$n_iter,
        backend = x$backend,
        memory_mb = x$memory_mb
      )
    }))

    memory_comparison <- memory_df %>%
      group_by(data_size, n_iter) %>%
      summarise(
        r_memory = memory_mb[backend == "r"],
        cpp_memory = memory_mb[backend == "cpp"],
        memory_ratio = memory_mb[backend == "r"] / memory_mb[backend == "cpp"],
        .groups = "drop"
      )

    cat("\n\nMemory Usage Summary:\n")
    cat("--------------------\n")
    print(memory_comparison)
  }

  # Process quality results
  if (length(results$quality) > 0) {
    quality_df <- do.call(rbind, lapply(results$quality, function(x) {
      data.frame(
        data_size = x$data_size,
        n_iter = x$n_iter,
        r_clusters = x$r_clusters,
        cpp_clusters = x$cpp_clusters,
        cluster_diff = abs(x$r_clusters - x$cpp_clusters),
        alpha_diff = abs(x$r_alpha - x$cpp_alpha),
        loglik_diff = abs(x$r_loglik - x$cpp_loglik)
      )
    }))

    cat("\n\nQuality Comparison Summary:\n")
    cat("--------------------------\n")
    print(quality_df)
    cat(sprintf("Average cluster difference: %.2f\n", mean(quality_df$cluster_diff)))
    cat(sprintf("Average alpha difference: %.4f\n", mean(quality_df$alpha_diff)))
  }

  return(list(speed = speed_df, memory = memory_comparison, quality = quality_df))
}

# Visualization Functions
create_benchmark_plots <- function(analyzed_results) {
  plots <- list()

  # Speed comparison plot
  if ("speed" %in% names(analyzed_results) && nrow(analyzed_results$speed) > 0) {
    speed_plot <- ggplot(analyzed_results$speed) +
      geom_point(aes(x = data_size, y = speedup, color = factor(n_iter)), size = 3) +
      geom_line(aes(x = data_size, y = speedup, color = factor(n_iter))) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_y_log10() +
      labs(
        title = "C++ vs R Speed Comparison",
        subtitle = "Higher values indicate C++ is faster",
        x = "Data Size",
        y = "Speedup Factor (log scale)",
        color = "Iterations"
      ) +
      theme_minimal()

    plots$speed <- speed_plot
  }

  # Memory comparison plot
  if ("memory" %in% names(analyzed_results) && nrow(analyzed_results$memory) > 0) {
    memory_plot <- ggplot(analyzed_results$memory) +
      geom_point(aes(x = data_size, y = memory_ratio, color = factor(n_iter)), size = 3) +
      geom_line(aes(x = data_size, y = memory_ratio, color = factor(n_iter))) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Memory Usage Comparison (R / C++)",
        subtitle = "Values > 1 indicate R uses more memory",
        x = "Data Size",
        y = "Memory Ratio",
        color = "Iterations"
      ) +
      theme_minimal()

    plots$memory <- memory_plot
  }

  return(plots)
}

# Quick Benchmark Function
quick_benchmark <- function(data_size = 200, n_iter = 100, times = 3) {
  cat("Quick Benchmark\n")
  cat("===============\n")
  cat(sprintf("Data size: %d, Iterations: %d, Replicates: %d\n\n", data_size, n_iter, times))

  data <- create_benchmark_data(data_size, 42)

  # Speed test
  cat("Running speed benchmark...\n")
  mb <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      dp <- DirichletProcessGaussian(data)
      Fit(dp, n_iter, progressBar = FALSE)
    },
    CPP = {
      set_use_cpp(TRUE)
      dp <- DirichletProcessGaussian(data)
      Fit(dp, n_iter, progressBar = FALSE)
    },
    times = times
  )

  print(mb)

  speedup <- median(mb$time[mb$expr == "R"]) / median(mb$time[mb$expr == "CPP"])
  cat(sprintf("\nSpeedup: %.2fx\n", speedup))

  # Quality comparison
  cat("\nRunning quality comparison...\n")
  quality <- compare_quality(data_size, n_iter)
  cat(sprintf("R clusters: %d, C++ clusters: %d\n", quality$r_clusters, quality$cpp_clusters))
  cat(sprintf("R alpha: %.3f, C++ alpha: %.3f\n", quality$r_alpha, quality$cpp_alpha))

  return(list(benchmark = mb, quality = quality, speedup = speedup))
}

# Example usage
if (interactive()) {
  cat("Dirichlet Process R vs C++ Benchmark\n")
  cat("====================================\n\n")

  cat("Available functions:\n")
  cat("- quick_benchmark(): Fast comparison with default parameters\n")
  cat("- run_comprehensive_benchmark(): Full benchmark suite\n")
  cat("- analyze_benchmark_results(): Process benchmark results\n")
  cat("- create_benchmark_plots(): Generate visualization plots\n\n")

  cat("Example usage:\n")
  cat("results <- quick_benchmark()\n")
  cat("# or\n")
  cat("full_results <- run_comprehensive_benchmark()\n")
  cat("analyzed <- analyze_benchmark_results(full_results)\n")
  cat("plots <- create_benchmark_plots(analyzed)\n")
}
