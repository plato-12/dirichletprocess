#' Fair component benchmarking for Weibull
#' @keywords internal
benchmark_weibull_components_fair <- function(n_obs = 1000, n_clusters = 3) {

  cat("Running fair component benchmarks for Weibull...\n")

  # Generate test data
  set.seed(123)
  shape_vec <- seq(1.5, 4.5, length.out = n_clusters)
  scale_vec <- seq(0.5, 2.0, length.out = n_clusters)
  data_list <- lapply(1:n_clusters, function(i) {
    rweibull(n_obs/n_clusters, shape = shape_vec[i], scale = scale_vec[i])
  })
  data_vec <- unlist(data_list)
  data <- matrix(data_vec, ncol = 1)

  # Initialize DP objects
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessWeibull(data)
  dp_r <- Initialise(dp_r, numInitialClusters = n_clusters)

  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessWeibull(data)
  dp_cpp <- Initialise(dp_cpp, numInitialClusters = n_clusters)

  # Fair benchmark: Run multiple iterations to amortize overhead
  n_iterations <- 100

  # 1. Likelihood calculation benchmark (batch mode)
  cat("  - Likelihood calculation (batch)...\n")

  # R version
  time_lik_r <- system.time({
    set_use_cpp(FALSE)
    for (iter in 1:n_iterations) {
      for (i in 1:nrow(data)) {
        Likelihood(dp_r$mixingDistribution,
                   data[i, , drop = FALSE],
                   dp_r$clusterParameters)
      }
    }
  })[["elapsed"]]

  # C++ version
  time_lik_cpp <- system.time({
    set_use_cpp(TRUE)
    for (iter in 1:n_iterations) {
      for (i in 1:nrow(data)) {
        Likelihood(dp_cpp$mixingDistribution,
                   data[i, , drop = FALSE],
                   dp_cpp$clusterParameters)
      }
    }
  })[["elapsed"]]

  # 2. Full MCMC iteration benchmark
  cat("  - Full MCMC iteration...\n")

  # R version
  time_iter_r <- system.time({
    set_use_cpp(FALSE)
    dp_r_temp <- dp_r
    for (iter in 1:10) {
      ClusterComponentUpdate(dp_r_temp)
      ClusterParameterUpdate(dp_r_temp)
      UpdateAlpha(dp_r_temp)
    }
  })[["elapsed"]]

  # C++ version
  time_iter_cpp <- system.time({
    set_use_cpp(TRUE)
    dp_cpp_temp <- dp_cpp
    for (iter in 1:10) {
      ClusterComponentUpdate(dp_cpp_temp)
      ClusterParameterUpdate(dp_cpp_temp)
      UpdateAlpha(dp_cpp_temp)
    }
  })[["elapsed"]]

  results <- data.frame(
    operation = c("Likelihood (per call)", "Full MCMC iteration"),
    R_time_ms = c(time_lik_r * 1000 / (n_iterations * n_obs),
                  time_iter_r * 1000 / 10),
    Cpp_time_ms = c(time_lik_cpp * 1000 / (n_iterations * n_obs),
                    time_iter_cpp * 1000 / 10),
    speedup = c(time_lik_r / time_lik_cpp,
                time_iter_r / time_iter_cpp)
  )

  return(results)
}
