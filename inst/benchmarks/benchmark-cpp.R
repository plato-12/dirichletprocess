# tests/benchmarks/benchmark-cpp.R

library(microbenchmark)
library(ggplot2)

benchmark_cpp_vs_r <- function(n_obs = c(100, 500, 1000),
                               n_iter = 1000) {
  results <- list()

  for (n in n_obs) {
    # Generate data
    data <- c(rnorm(n/2, -2, 1), rnorm(n/2, 2, 1))

    # Create DP objects
    dp_r <- DirichletProcessGaussian(data)
    dp_cpp <- DirichletProcessGaussian(data)

    # Benchmark
    mb <- microbenchmark(
      R = {
        set_use_cpp(FALSE)
        Fit(dp_r, n_iter, n_iter/2, 1)
      },
      CPP = {
        set_use_cpp(TRUE)
        Fit(dp_cpp, n_iter, n_iter/2, 1)
      },
      times = 10
    )

    results[[paste0("n", n)]] <- mb
  }

  # Plot results
  plot_data <- do.call(rbind, lapply(names(results), function(name) {
    df <- as.data.frame(results[[name]])
    df$n_obs <- as.numeric(gsub("n", "", name))
    df
  }))

  ggplot(plot_data, aes(x = as.factor(n_obs), y = time/1e9, fill = expr)) +
    geom_boxplot() +
    labs(x = "Number of Observations", y = "Time (seconds)",
         title = "R vs C++ Implementation Performance") +
    scale_fill_manual(values = c("R" = "coral", "CPP" = "steelblue"))
}

# Run benchmark
benchmark_cpp_vs_r()
