# inst/benchmarks/comprehensive_benchmark.R
library(dirichletprocess)
library(microbenchmark)
library(ggplot2)

run_comprehensive_benchmark <- function(output_dir = "benchmark_results") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Configuration
  sample_sizes <- c(100, 500, 1000, 5000)
  distributions <- c("gaussian", "beta", "mvnormal")
  iterations <- 10  # Number of MCMC iterations for Fit benchmarks

  # Results data frame
  results <- data.frame()

  # 1. Likelihood Function Benchmarks
  cat("Running Likelihood benchmarks...\n")

  for (dist in distributions) {
    for (n in sample_sizes) {
      # Setup test data and objects
      if (dist == "gaussian") {
        x <- rnorm(n)
        mdObj <- MixingDistribution("normal", c(0,1,1,1), "conjugate")
        theta <- list(array(0, dim=c(1,1,1)), array(1, dim=c(1,1,1)))
      } else if (dist == "beta") {
        x <- rbeta(n, 2, 5)
        mdObj <- BetaMixtureCreate(c(2,8), c(1,1), 1)
        theta <- list(array(0.5, dim=c(1,1,1)), array(2, dim=c(1,1,1)))
      } else if (dist == "mvnormal") {
        x <- matrix(rnorm(n*2), ncol=2)
        priorParams <- list(mu0=c(0,0), Lambda=diag(2), kappa0=1, nu=2)
        mdObj <- MvnormalCreate(priorParams)
        theta <- list(
          array(c(0,0), dim=c(1,2,1)),
          array(diag(2), dim=c(2,2,1))
        )
      }

      # R implementation
      set_use_cpp(FALSE)
      r_bench <- microbenchmark(
        Likelihood(mdObj, x, theta),
        times = 20
      )

      # C++ implementation
      set_use_cpp(TRUE)
      cpp_bench <- microbenchmark(
        Likelihood(mdObj, x, theta),
        times = 20
      )

      # Record results
      results <- rbind(results, data.frame(
        function_name = "Likelihood",
        distribution = dist,
        n = n,
        implementation = "R",
        median_time_ms = median(r_bench$time) / 1e6,
        mean_time_ms = mean(r_bench$time) / 1e6
      ))

      results <- rbind(results, data.frame(
        function_name = "Likelihood",
        distribution = dist,
        n = n,
        implementation = "C++",
        median_time_ms = median(cpp_bench$time) / 1e6,
        mean_time_ms = mean(cpp_bench$time) / 1e6
      ))
    }
  }

  # 2. ClusterComponentUpdate Benchmarks
  cat("Running ClusterComponentUpdate benchmarks...\n")

  for (dist in distributions) {
    for (n in sample_sizes[1:2]) {  # Use smaller sample sizes
      # Create and initialize DP object
      if (dist == "gaussian") {
        dp <- DirichletProcessGaussian(rnorm(n))
      } else if (dist == "beta") {
        dp <- DirichletProcessBeta(rbeta(n, 2, 5), 1, verbose=FALSE)
      } else if (dist == "mvnormal") {
        dp <- DirichletProcessMvnormal(matrix(rnorm(n*2), ncol=2))
      }

      # Run a few iterations to establish realistic clustering
      dp <- Fit(dp, 5, progressBar=FALSE)

      # R implementation
      set_use_cpp(FALSE)
      r_bench <- microbenchmark(
        ClusterComponentUpdate(dp),
        times = 5
      )

      # C++ implementation (when available)
      set_use_cpp(TRUE)
      tryCatch({
        cpp_bench <- microbenchmark(
          ClusterComponentUpdate(dp),
          times = 5
        )

        results <- rbind(results, data.frame(
          function_name = "ClusterComponentUpdate",
          distribution = dist,
          n = n,
          implementation = "C++",
          median_time_ms = median(cpp_bench$time) / 1e6,
          mean_time_ms = mean(cpp_bench$time) / 1e6
        ))
      }, error = function(e) {
        cat("C++ implementation not available for", dist, "ClusterComponentUpdate\n")
      })

      # Always record R results
      results <- rbind(results, data.frame(
        function_name = "ClusterComponentUpdate",
        distribution = dist,
        n = n,
        implementation = "R",
        median_time_ms = median(r_bench$time) / 1e6,
        mean_time_ms = mean(r_bench$time) / 1e6
      ))
    }
  }

  # 3. ClusterParameterUpdate Benchmarks
  cat("Running ClusterParameterUpdate benchmarks...\n")

  # (Similar structure to ClusterComponentUpdate benchmarks)

  # 4. Fit Benchmarks (for small dataset only)
  cat("Running Fit benchmarks...\n")

  for (dist in distributions) {
    n <- sample_sizes[1]  # Use smallest sample size

    # Create DP object
    if (dist == "gaussian") {
      dp <- DirichletProcessGaussian(rnorm(n))
    } else if (dist == "beta") {
      dp <- DirichletProcessBeta(rbeta(n, 2, 5), 1, verbose=FALSE)
    } else if (dist == "mvnormal") {
      dp <- DirichletProcessMvnormal(matrix(rnorm(n*2), ncol=2))
    }

    # R implementation
    set_use_cpp(FALSE)
    r_time <- system.time(Fit(dp, iterations, progressBar=FALSE))

    # C++ implementation (when available)
    set_use_cpp(TRUE)
    tryCatch({
      cpp_time <- system.time(Fit(dp, iterations, progressBar=FALSE))

      results <- rbind(results, data.frame(
        function_name = "Fit",
        distribution = dist,
        n = n,
        implementation = "C++",
        median_time_ms = cpp_time["elapsed"] * 1000 / iterations,
        mean_time_ms = cpp_time["elapsed"] * 1000 / iterations
      ))
    }, error = function(e) {
      cat("C++ implementation not available for", dist, "Fit\n")
    })

    # Always record R results
    results <- rbind(results, data.frame(
      function_name = "Fit",
      distribution = dist,
      n = n,
      implementation = "R",
      median_time_ms = r_time["elapsed"] * 1000 / iterations,
      mean_time_ms = r_time["elapsed"] * 1000 / iterations
    ))
  }

  # Save results
  write.csv(results, file.path(output_dir, "benchmark_results.csv"), row.names = FALSE)

  # Generate summary plots
  generate_benchmark_plots(results, output_dir)

  # Reset to default
  set_use_cpp(FALSE)

  return(results)
}

# Helper function to create plots
generate_benchmark_plots <- function(results, output_dir) {
  # Likelihood plots
  likelihood_results <- subset(results, function_name == "Likelihood")

  p1 <- ggplot(likelihood_results,
               aes(x = n, y = median_time_ms, color = implementation)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ distribution, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Likelihood Performance by Distribution",
      x = "Sample Size (log scale)",
      y = "Median Time (ms, log scale)",
      color = "Implementation"
    ) +
    theme_minimal()

  ggsave(file.path(output_dir, "likelihood_benchmark.png"), p1, width = 10, height = 6)

  # Implementation comparison for each function
  for (func in unique(results$function_name)) {
    func_results <- subset(results, function_name == func)

    if (nrow(func_results) > 0) {
      p <- ggplot(func_results,
                  aes(x = n, y = median_time_ms, color = implementation,
                      shape = distribution)) +
        geom_line() +
        geom_point(size = 3) +
        scale_x_log10() +
        scale_y_log10() +
        labs(
          title = paste(func, "Performance"),
          x = "Sample Size (log scale)",
          y = "Median Time (ms, log scale)",
          color = "Implementation",
          shape = "Distribution"
        ) +
        theme_minimal()

      ggsave(file.path(output_dir, paste0(func, "_benchmark.png")), p, width = 10, height = 6)
    }
  }

  # Speedup ratios
  calculate_speedups <- function(df) {
    # Reshape to wide format
    wide <- reshape(df,
                    idvar = c("function_name", "distribution", "n"),
                    timevar = "implementation",
                    direction = "wide")

    if (ncol(wide) >= 4) {  # Make sure we have both R and C++ results
      wide$speedup <- wide$median_time_ms.R / wide$median_time_ms.C++
        return(wide)
    }
    return(NULL)
  }

  # Calculate speedups by function and distribution
  speedups <- NULL
  for (func in unique(results$function_name)) {
    for (dist in unique(results$distribution)) {
      subset_data <- results[results$function_name == func &
                               results$distribution == dist, ]

      if (nrow(subset_data) > 0) {
        speedup_data <- calculate_speedups(subset_data)
        if (!is.null(speedup_data)) {
          speedups <- rbind(speedups, speedup_data)
        }
      }
    }
  }

  if (!is.null(speedups) && nrow(speedups) > 0) {
    p <- ggplot(speedups, aes(x = n, y = speedup, color = distribution)) +
      geom_line() +
      geom_point() +
      facet_wrap(~ function_name, scales = "free_y") +
      scale_x_log10() +
      labs(
        title = "C++ Speedup Ratio by Function",
        x = "Sample Size (log scale)",
        y = "Speedup Ratio (R/C++)",
        color = "Distribution"
      ) +
      theme_minimal()

    ggsave(file.path(output_dir, "speedup_ratios.png"), p, width = 10, height = 6)
  }
}

# Run the benchmarks
results <- run_comprehensive_benchmark()
