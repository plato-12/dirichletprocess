# Define a log file path with a timestamp for uniqueness
log_file <- file.path("analysis/benchmark_results",
                      paste0("benchmark_log_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt"))

# Redirect all console output to the log file
sink(log_file)

# --- Your Benchmark Code ---
# It's good practice to print session info for reproducibility
print(sessionInfo())
cat("\n\n")

# Now run the benchmark
library(dirichletprocess)

# R vs C++ Benchmark Function
run_benchmark <- function() {
  cat("=== R vs C++ Performance Benchmark ===\n")

  sizes <- c(50, 100, 200)
  iterations <- 10

  for (n in sizes) {
    cat(sprintf("\nBenchmarking with %d data points, %d iterations...\n", n, iterations))

    # Create test data with clear clusters
    set.seed(123)
    y <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 2, 0.5))

    # R implementation benchmark
    dp_r <- DirichletProcessGaussian(y)
    time_r <- system.time({
      for (i in 1:iterations) {
        dp_r <- ClusterComponentUpdate(dp_r)
        dp_r <- ClusterParameterUpdate(dp_r)
      }
    })

    # C++ implementation benchmark
    dp_cpp <- DirichletProcessGaussian(y)
    dp_cpp$clusterLabels <- dp_cpp$clusterLabels - 1  # Convert to 0-indexed

    time_cpp <- system.time({
      for (i in 1:iterations) {
        # C++ cluster component update
        update_result <- conjugate_cluster_component_update_cpp(dp_cpp)
        dp_cpp$clusterLabels <- update_result$clusterLabels
        dp_cpp$pointsPerCluster <- update_result$pointsPerCluster
        dp_cpp$numberClusters <- update_result$numberClusters
        dp_cpp$clusterParameters <- update_result$clusterParameters

        # C++ cluster parameter update
        dp_cpp$clusterParameters <- conjugate_cluster_parameter_update_cpp(dp_cpp)
      }
    })

    # Convert back to 1-indexed for comparison
    dp_cpp$clusterLabels <- dp_cpp$clusterLabels + 1

    # Results
    speedup <- time_r["elapsed"] / time_cpp["elapsed"]
    cat(sprintf("  R implementation:   %.3f seconds\n", time_r["elapsed"]))
    cat(sprintf("  C++ implementation: %.3f seconds\n", time_cpp["elapsed"]))
    cat(sprintf("  Speedup: %.1fx\n", speedup))
    cat(sprintf("  R clusters: %d, C++ clusters: %d\n",
                dp_r$numberClusters, dp_cpp$numberClusters))
  }
}

# Run the benchmark
run_benchmark()

# --- End of Benchmark Code ---

# Stop redirecting output and close the file connection
sink()

# Optional: Print a message to the console to confirm completion
message("Benchmark complete. Log saved to: ", log_file)
