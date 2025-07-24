# benchmark/manual_mcmc_performance.R

benchmark_manual_mcmc <- function() {
  cat("\nBenchmarking manual MCMC interface...\n")

  results <- list()

  for (dist in c("normal", "exponential", "mvnormal")) {
    test_data <- generate_test_data(dist, n = 1000)
    dp <- create_dp_object(dist, test_data)

    # Benchmark individual step functions
    runner <- CppMCMCRunner$new(dp)

    step_times <- microbenchmark(
      assignments = runner$step_assignments(),
      parameters = runner$step_parameters(),
      concentration = runner$step_concentration(),
      times = 100
    )

    # Compare with full step
    full_step_time <- microbenchmark(
      full = {
        runner$step_assignments()
        runner$step_parameters()
        runner$step_concentration()
      },
      times = 100
    )

    results[[dist]] <- list(
      step_times = step_times,
      full_step_time = full_step_time,
      overhead = median(full_step_time$time) -
        sum(sapply(split(step_times$time, step_times$expr), median))
    )
  }

  return(results)
}

# Compare manual MCMC vs Fit() performance
benchmark_manual_vs_fit <- function() {
  cat("\nComparing manual MCMC vs Fit() performance...\n")

  distributions <- c("normal", "exponential", "beta", "weibull")
  results <- list()

  for (dist in distributions) {
    cat("  Testing", dist, "...")

    test_data <- generate_test_data(dist, n = 500)
    iterations <- 100

    # Benchmark Fit() approach
    fit_time <- microbenchmark(
      fit = {
        dp <- create_dp_object(dist, test_data)
        Fit(dp, its = iterations)
      },
      times = 10
    )

    # Benchmark manual approach
    manual_time <- microbenchmark(
      manual = {
        dp <- create_dp_object(dist, test_data)
        runner <- CppMCMCRunner$new(dp)
        for (i in 1:iterations) {
          runner$step_assignments()
          runner$step_parameters()
          runner$step_concentration()
        }
      },
      times = 10
    )

    overhead_pct <- (median(manual_time$time) - median(fit_time$time)) /
      median(fit_time$time) * 100

    results[[dist]] <- list(
      fit_time = median(fit_time$time),
      manual_time = median(manual_time$time),
      overhead_percent = overhead_pct
    )

    cat(" Overhead:", round(overhead_pct, 1), "%\n")
  }

  return(results)
}

# Benchmark advanced features
benchmark_advanced_features <- function() {
  cat("\nBenchmarking advanced CppMCMCRunner features...\n")

  test_data <- generate_test_data("normal", n = 1000)
  dp <- DirichletProcessGaussian(test_data)
  runner <- CppMCMCRunner$new(dp)

  # Warm up
  for (i in 1:10) {
    runner$step_assignments()
    runner$step_parameters()
  }

  # Benchmark various operations
  operations <- microbenchmark(
    get_state = runner$get_state(),
    get_n_clusters = runner$get_n_clusters(),
    sample_predictive_10 = runner$sample_predictive(n = 10),
    sample_predictive_100 = runner$sample_predictive(n = 100),
    set_temperature = runner$set_temperature(0.5),
    get_temperature = runner$get_temperature(),
    times = 100
  )

  print(operations)

  return(operations)
}

# Benchmark temperature-controlled MCMC
benchmark_temperature_mcmc <- function() {
  cat("\nBenchmarking temperature-controlled MCMC...\n")

  test_data <- generate_test_data("exponential", n = 500)
  dp <- DirichletProcessExponential(test_data)

  temperatures <- c(1.0, 0.8, 0.5, 0.2)
  results <- list()

  for (temp in temperatures) {
    runner <- CppMCMCRunner$new(dp)
    runner$set_temperature(temp)

    time_result <- microbenchmark(
      tempered_step = {
        runner$step_assignments()
        runner$step_parameters()
        runner$step_concentration()
      },
      times = 50
    )

    results[[as.character(temp)]] <- median(time_result$time)
  }

  # Plot temperature vs time
  temp_df <- data.frame(
    temperature = temperatures,
    time = unlist(results)
  )

  p <- ggplot(temp_df, aes(x = temperature, y = time)) +
    geom_line() +
    geom_point() +
    labs(title = "MCMC Step Time vs Temperature",
         x = "Temperature", y = "Time (nanoseconds)") +
    theme_minimal()

  print(p)

  return(results)
}

# Profile memory usage of manual MCMC
profile_manual_mcmc_memory <- function() {
  if (!requireNamespace("profmem", quietly = TRUE)) {
    warning("profmem package not installed")
    return(NULL)
  }

  cat("\nProfiling manual MCMC memory usage...\n")

  test_data <- generate_test_data("mvnormal", n = 1000)
  dp <- DirichletProcessMvnormal(test_data)

  # Profile creation
  creation_mem <- profmem::profmem({
    runner <- CppMCMCRunner$new(dp)
  })

  runner <- CppMCMCRunner$new(dp)

  # Profile single iteration
  iteration_mem <- profmem::profmem({
    runner$step_assignments()
    runner$step_parameters()
    runner$step_concentration()
  })

  # Profile state extraction
  state_mem <- profmem::profmem({
    state <- runner$get_state()
  })

  # Profile predictive sampling
  predictive_mem <- profmem::profmem({
    samples <- runner$sample_predictive(n = 100)
  })

  results <- list(
    creation = profmem::total(creation_mem),
    iteration = profmem::total(iteration_mem),
    state_extraction = profmem::total(state_mem),
    predictive_sampling = profmem::total(predictive_mem)
  )

  cat("Memory usage (bytes):\n")
  print(results)

  return(results)
}
