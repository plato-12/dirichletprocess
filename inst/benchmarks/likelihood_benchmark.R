# Benchmark script for likelihood calculations
library(dirichletprocess)
library(microbenchmark)

# Set up test data
set.seed(123)
sample_sizes <- c(10, 100, 1000, 10000)

# Run benchmarks
results <- data.frame()

for (n in sample_sizes) {
  x <- rnorm(n)
  mdObj <- MixingDistribution("normal", c(0,1,1,1), "conjugate")
  theta <- list(array(0, dim=c(1,1,1)), array(1, dim=c(1,1,1)))

  # R implementation
  set_use_cpp(FALSE)
  r_bench <- microbenchmark(
    Likelihood(mdObj, x, theta),
    times = 10
  )

  # C++ implementation
  set_use_cpp(TRUE)
  cpp_bench <- microbenchmark(
    Likelihood(mdObj, x, theta),
    times = 10
  )

  # Add results
  results <- rbind(results, data.frame(
    n = n,
    implementation = "R",
    median_time_ms = median(r_bench$time) / 1e6
  ))

  results <- rbind(results, data.frame(
    n = n,
    implementation = "C++",
    median_time_ms = median(cpp_bench$time) / 1e6
  ))
}

# Print results
print(results)

# Reset to default
set_use_cpp(FALSE)

# Plot results if ggplot2 is available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  ggplot(results, aes(x = n, y = median_time_ms, color = implementation)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Likelihood Calculation Performance",
      x = "Sample Size (log scale)",
      y = "Median Time (ms, log scale)",
      color = "Implementation"
    ) +
    theme_minimal()

  ggsave("likelihood_benchmark.png", width = 8, height = 6)
}
