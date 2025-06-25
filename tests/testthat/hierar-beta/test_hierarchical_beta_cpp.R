context("Hierarchical Beta C++ Implementation Tests")

test_that("Hierarchical Beta MCMC runs without errors", {
  skip_if_not(exists("_dirichletprocess_run_hierarchical_mcmc_cpp"))

  # Generate hierarchical data
  set.seed(123)
  true_mu <- c(0.3, 0.7)
  true_tau <- c(20, 15)

  # Two datasets with shared clusters
  data1 <- c(
    rbeta(50, true_mu[1] * true_tau[1], (1 - true_mu[1]) * true_tau[1]),
    rbeta(50, true_mu[2] * true_tau[2], (1 - true_mu[2]) * true_tau[2])
  )

  data2 <- c(
    rbeta(60, true_mu[1] * true_tau[1], (1 - true_mu[1]) * true_tau[1]),
    rbeta(40, true_mu[2] * true_tau[2], (1 - true_mu[2]) * true_tau[2])
  )

  # Create hierarchical DP
  dp_list <- DirichletProcessHierarchicalBeta(
    list(data1, data2),
    maxY = 1,
    hyperPriorParameters = c(1, 0.01),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4)
  )

  # Run with C++ implementation
  if (can_use_hierarchical_cpp(dp_list)) {
    result <- run_hierarchical_mcmc_cpp(
      dp_list,
      n_iter = 100,
      n_burn = 20,
      progress_bar = FALSE
    )

    expect_true(is.list(result))
    expect_true(length(result$indDP) == 2)
    expect_true(length(result$gammaValues) > 0)
    expect_true(all(result$gammaValues > 0))

    # Check that clusters were found
    for (dp in result$indDP) {
      expect_true(dp$numberClusters >= 1)
      expect_true(dp$numberClusters <= 10)
    }
  }
})

test_that("Hierarchical Beta parameters are updated correctly", {
  # Create dummy data first to check if we can use hierarchical CPP
  dummy_data <- list(rbeta(10, 2, 8))
  dummy_dp <- DirichletProcessHierarchicalBeta(
    dummy_data,
    maxY = 1,
    hyperPriorParameters = c(1, 0.01)
  )

  skip_if_not(can_use_hierarchical_cpp(dummy_dp))

  # Create simple test case
  data_list <- list(
    rbeta(50, 2, 8),
    rbeta(50, 2, 8)
  )

  dp_list <- DirichletProcessHierarchicalBeta(
    data_list,
    maxY = 1,
    hyperPriorParameters = c(1, 0.01)
  )

  # Initial values
  initial_gamma <- dp_list$gamma

  # Run MCMC
  result <- run_hierarchical_mcmc_cpp(dp_list, n_iter = 50, progress_bar = FALSE)

  # Check updates occurred
  expect_true(result$gamma != initial_gamma)
  expect_true(length(result$globalParameters) > 0)
})

test_that("Performance improvement over R implementation", {
  skip_if_not_installed("microbenchmark")

  # Create dummy data first to check if we can use hierarchical CPP
  dummy_data <- list(rbeta(10, 2, 5))
  dummy_dp <- DirichletProcessHierarchicalBeta(
    dummy_data,
    maxY = 1
  )

  skip_if_not(can_use_hierarchical_cpp(dummy_dp))

  # Generate test data
  set.seed(42)
  data_list <- lapply(1:3, function(i) rbeta(100, 2, 5))

  dp_list <- DirichletProcessHierarchicalBeta(
    data_list,
    maxY = 1
  )

  # Benchmark
  mb_result <- microbenchmark::microbenchmark(
    R_Backend = {
      set_use_cpp(FALSE)
      Fit(dp_list, 20, progressBar = FALSE)
    },
    CPP_Backend = {
      result <- run_hierarchical_mcmc_cpp(dp_list, n_iter = 20, progress_bar = FALSE)
    },
    times = 3
  )

  # C++ should be faster
  r_time <- median(mb_result$time[mb_result$expr == "R_Backend"])
  cpp_time <- median(mb_result$time[mb_result$expr == "CPP_Backend"])
  speedup <- r_time / cpp_time

  expect_true(speedup > 5, info = sprintf("Speedup: %.1fx", speedup))
})
