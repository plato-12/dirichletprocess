# Save as profiling/benchmarks.R
library(dirichletprocess)
library(microbenchmark)

benchmark_gaussian <- function(n_obs = 500, n_iter = 10) {
  set.seed(123)
  y <- rnorm(n_obs)

  # Create DP object first to focus on the fitting time
  dp <- DirichletProcessGaussian(y)

  # Benchmark just the fitting part
  results <- microbenchmark(
    Fit(dp, n_iter, progressBar = FALSE),
    times = 5
  )

  print(results)
  return(results)
}

benchmark_gaussian()
