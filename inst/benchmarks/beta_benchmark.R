library(dirichletprocess)
library(microbenchmark)

# Beta distribution benchmark
run_beta_benchmark <- function(n_data = c(50, 100, 200),
                               n_iterations = 20) {

  results <- data.frame()

  for (n in n_data) {
    cat("Benchmarking Beta with", n, "data points...\n")

    # Generate test data
    set.seed(123)
    y <- rbeta(n, 3, 7)

    # R implementation
    enable_cpp_samplers(FALSE)
    dp_r <- DirichletProcessBeta(y, 1, verbose = FALSE)

    time_r <- microbenchmark(
      Fit(dp_r, n_iterations, progressBar = FALSE),
      times = 5
    )

    # C++ implementation
    enable_cpp_samplers(TRUE)
    dp_cpp <- DirichletProcessBeta(y, 1, verbose = FALSE)

    time_cpp <- microbenchmark(
      Fit(dp_cpp, n_iterations, progressBar = FALSE),
      times = 5
    )

    results <- rbind(results, data.frame(
      n = n,
      implementation = c("R", "C++"),
      median_time = c(median(time_r$time)/1e9,
                      median(time_cpp$time)/1e9),
      speedup = median(time_r$time) / median(time_cpp$time)
    ))
  }

  return(results)
}

# Run benchmark
results <- run_beta_benchmark()
print(results)
