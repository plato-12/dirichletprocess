# benchmarks/benchmark_mvnormal.R
library(dirichletprocess)
library(mvtnorm)
library(microbenchmark)

benchmark_mvnormal <- function(n_obs = 100, d = 2, n_iter = 100) {
  # Generate test data
  set.seed(123)
  y <- mvtnorm::rmvnorm(n_obs, rep(0, d), diag(d))

  # Benchmark
  mb <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      dp <- DirichletProcessMvnormal(y)
      Fit(dp, n_iter, progressBar = FALSE)
    },
    Cpp = {
      set_use_cpp(TRUE)
      dp <- DirichletProcessMvnormal(y)
      Fit(dp, n_iter, progressBar = FALSE)
    },
    times = 5
  )

  print(mb)

  # Calculate speedup
  times <- summary(mb)
  speedup <- times$mean[1] / times$mean[2]
  cat(sprintf("\nSpeedup: %.1fx\n", speedup))

  return(mb)
}

# Run benchmark
benchmark_mvnormal(n_obs = 200, d = 3, n_iter = 100)
