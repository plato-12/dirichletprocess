# tests/testthat/test-cpp-consistency.R

library(testthat)
library(dirichletprocess)

# Tolerance levels for statistical tests
ALPHA_TOLERANCE <- 0.05      # Mean alpha difference
CLUSTER_TOLERANCE <- 0.1     # Mean cluster count difference
LIKELIHOOD_CORR_MIN <- 0.95  # Minimum likelihood correlation
PARAM_TOLERANCE <- 0.05      # Parameter estimate differences

# Main consistency validation function
validate_r_cpp_consistency <- function(distribution_type,
                                       test_data,
                                       iterations = 100,
                                       n_runs = 5,
                                       seed = 12345) {

  consistency_results <- list()

  for (run in 1:n_runs) {
    current_seed <- seed + run - 1

    # R implementation
    set.seed(current_seed)
    set_use_cpp(FALSE)
    dp_r <- create_dp_object(distribution_type, test_data)
    dp_r <- Fit(dp_r, its = iterations)

    # C++ implementation
    set.seed(current_seed)
    set_use_cpp(TRUE)
    dp_cpp <- create_dp_object(distribution_type, test_data)
    dp_cpp <- Fit(dp_cpp, its = iterations)

    # Extract statistics
    r_stats <- extract_dp_statistics(dp_r)
    cpp_stats <- extract_dp_statistics(dp_cpp)

    # Statistical consistency checks
    consistency_results[[run]] <- list(
      alpha_mean_diff = abs(r_stats$alpha_mean - cpp_stats$alpha_mean),
      alpha_sd_diff = abs(r_stats$alpha_sd - cpp_stats$alpha_sd),
      cluster_count_diff = abs(r_stats$mean_clusters - cpp_stats$mean_clusters),
      likelihood_correlation = cor(dp_r$likelihoodChain, dp_cpp$likelihoodChain),
      param_max_diff = max(abs(unlist(r_stats$param_means) - unlist(cpp_stats$param_means))),
      runtime_r = r_stats$runtime,
      runtime_cpp = cpp_stats$runtime
    )
  }

  # Aggregate results
  aggregate_consistency_results(consistency_results)
}

# Helper function to extract statistics
extract_dp_statistics <- function(dp_obj) {
  start_time <- Sys.time()

  list(
    alpha_mean = mean(dp_obj$alphaChain),
    alpha_sd = sd(dp_obj$alphaChain),
    mean_clusters = mean(sapply(dp_obj$labelsChain, function(x) length(unique(x)))),
    param_means = lapply(dp_obj$clusterParametersChain, function(params) {
      if (is.list(params)) {
        lapply(params, mean)
      } else {
        mean(params)
      }
    }),
    runtime = as.numeric(Sys.time() - start_time)
  )
}

# Aggregate multiple runs
aggregate_consistency_results <- function(results) {
  list(
    alpha_mean_diff = mean(sapply(results, `[[`, "alpha_mean_diff")),
    alpha_sd_diff = mean(sapply(results, `[[`, "alpha_sd_diff")),
    cluster_count_diff = mean(sapply(results, `[[`, "cluster_count_diff")),
    likelihood_correlation = mean(sapply(results, `[[`, "likelihood_correlation")),
    param_max_diff = mean(sapply(results, `[[`, "param_max_diff")),
    speedup_factor = mean(sapply(results, function(x) x$runtime_r / x$runtime_cpp)),
    all_runs = results
  )
}

# Basic test to ensure framework is working
test_that("Consistency validation framework works", {
  test_data <- rnorm(50)
  results <- validate_r_cpp_consistency("normal", test_data, iterations = 10, n_runs = 2)

  expect_type(results, "list")
  expect_true("alpha_mean_diff" %in% names(results))
  expect_true("speedup_factor" %in% names(results))
  expect_true(results$speedup_factor > 0)
})
