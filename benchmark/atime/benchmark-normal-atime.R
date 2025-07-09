# benchmark_normal_atime.R
# Benchmark Normal Distribution using atime package
# Based on Neal (2000) and Escobar & West (1995) algorithms

library(dirichletprocess)
library(atime)
library(ggplot2)

# ==============================================================================
# Setup Functions
# ==============================================================================

#' Generate synthetic data for benchmarking
#' Creates mixture of normal distributions with known clusters
generate_test_data <- function(n, k_clusters = 3, seed = 42) {
  set.seed(seed)

  # Equal sized clusters
  cluster_sizes <- rep(n %/% k_clusters, k_clusters)
  cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)

  # Well-separated cluster means
  means <- seq(-3, 3, length.out = k_clusters)
  sds <- rep(0.8, k_clusters)

  data <- numeric(n)
  idx <- 1

  for (i in 1:k_clusters) {
    cluster_data <- rnorm(cluster_sizes[i], mean = means[i], sd = sds[i])
    data[idx:(idx + cluster_sizes[i] - 1)] <- cluster_data
    idx <- idx + cluster_sizes[i]
  }

  return(data)
}

# ==============================================================================
# Main atime Benchmark
# ==============================================================================

#' Run atime benchmark comparing R and C++ implementations
run_atime_benchmark <- function() {

  cat("========================================\n")
  cat("Normal Distribution DP Benchmark (atime)\n")
  cat("========================================\n\n")

  # Check C++ availability
  cpp_available <- exists("_dirichletprocess_run_mcmc_cpp") ||
    (exists("get_cpp_status") && get_cpp_status()$mcmc_runner)

  if (!cpp_available) {
    cat("⚠️  C++ implementation not available\n")
    cat("Only R implementation will be benchmarked\n\n")
  } else {
    cat("✓ C++ implementation available\n\n")
  }

  # Define expression list for atime
  expr_list <- list()

  # R implementation expression
  expr_list$R_implementation <- quote({
    set_use_cpp(FALSE)
    set.seed(123)
    dp <- DirichletProcessGaussian(data)
    dp <- Fit(dp, iterations, progressBar = FALSE)
  })

  # C++ implementation expression (if available)
  if (cpp_available) {
    expr_list$Cpp_implementation <- quote({
      set_use_cpp(TRUE)
      set.seed(123)
      dp <- DirichletProcessGaussian(data)
      dp <- Fit(dp, iterations, progressBar = FALSE)
    })
  }

  # Run atime benchmark
  atime_result <- atime::atime(
    N = as.integer(10^seq(1.5, 3.5, by = 0.25)),  # 30 to 3000+ observations
    setup = {
      data <- generate_test_data(N, k_clusters = 3)
      iterations <- 100  # Fixed number of iterations
    },
    expr.list = expr_list,
    seconds.limit = 1000000,  # Stop if any expression takes more than 10 seconds
    verbose = TRUE
  )

  # Print summary
  print(atime_result)

  # Create and save plot
  p <- plot(atime_result)
  print(p)

  # Save results
  cat("\nSaving results...\n")
  save(atime_result, file = "atime_normal_results.RData")
  ggsave("atime_normal_benchmark.png", plot = p, width = 10, height = 8, dpi = 300)

  return(atime_result)
}

# ==============================================================================
# Component-level Benchmarking with atime
# ==============================================================================

#' Benchmark individual components using atime
benchmark_components <- function() {

  cat("\n\n==========================================\n")
  cat("Component-level Benchmarking (atime)\n")
  cat("==========================================\n\n")

  # Check if component benchmarking is available
  if (!exists("benchmark_cpp_components")) {
    cat("Component benchmarking not available\n")
    return(NULL)
  }

  # Define components to benchmark
  components <- c("ClusterAssignment", "ParameterUpdate", "Likelihood")

  results <- list()

  for (comp in components) {
    cat(sprintf("\nBenchmarking %s...\n", comp))

    # R implementation
    expr_list_comp <- list()
    expr_list_comp[[paste0(comp, "_R")]] <- substitute({
      set_use_cpp(FALSE)
      dp <- DirichletProcessGaussian(data)
      # Initialize
      dp <- InitialiseClusters(dp)

      # Run specific component multiple times
      for (i in 1:10) {
        if (comp == "ClusterAssignment") {
          dp <- ClusterComponentUpdate(dp)
        } else if (comp == "ParameterUpdate") {
          dp <- ClusterParameterUpdate(dp)
        } else if (comp == "Likelihood") {
          ll <- sum(log(Likelihood(dp$mixingDistribution, data, dp$clusterParameters)))
        }
      }
    }, list(comp = comp))

    # C++ implementation (if available)
    if (exists("benchmark_cpp_components")) {
      expr_list_comp[[paste0(comp, "_Cpp")]] <- substitute({
        set_use_cpp(TRUE)
        dp <- DirichletProcessGaussian(data)
        result <- benchmark_cpp_components(dp, comp, times = 10)
      }, list(comp = comp))
    }

    # Run atime for this component
    comp_result <- atime::atime(
      N = as.integer(10^seq(2, 3.5, by = 0.5)),
      setup = {
        data <- generate_test_data(N)
      },
      expr.list = expr_list_comp,
      seconds.limit = 5
    )

    results[[comp]] <- comp_result
    print(comp_result)
  }

  return(results)
}

# ==============================================================================
# Memory Scaling Analysis
# ==============================================================================

#' Analyze memory usage scaling with atime
analyze_memory_scaling <- function() {

  cat("\n\n========================================\n")
  cat("Memory Scaling Analysis (atime)\n")
  cat("========================================\n\n")

  # Memory tracking expressions
  expr_list_mem <- list()

  # R implementation with memory tracking
  expr_list_mem$R_memory <- quote({
    set_use_cpp(FALSE)
    gc(reset = TRUE)

    dp <- DirichletProcessGaussian(data)
    dp <- Fit(dp, 50, progressBar = FALSE)

    mem_used <- gc()[2, 2]  # Total memory in MB
  })

  # C++ implementation with memory tracking
  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    expr_list_mem$Cpp_memory <- quote({
      set_use_cpp(TRUE)
      gc(reset = TRUE)

      if (exists("clear_memory_tracking")) clear_memory_tracking()

      dp <- DirichletProcessGaussian(data)
      dp <- Fit(dp, 50, progressBar = FALSE)

      mem_used <- gc()[2, 2]  # Total memory in MB
    })
  }

  # Run memory benchmark
  mem_result <- atime::atime(
    N = as.integer(10^seq(2, 3.5, by = 0.5)),
    setup = {
      data <- generate_test_data(N)
    },
    expr.list = expr_list_mem,
    seconds.limit = 10
  )

  print(mem_result)

  return(mem_result)
}

# ==============================================================================
# Comparison with Different Prior Settings
# ==============================================================================

#' Benchmark with different prior configurations
benchmark_prior_variations <- function() {

  cat("\n\n========================================\n")
  cat("Prior Configuration Benchmarks (atime)\n")
  cat("========================================\n\n")

  # Different prior configurations
  prior_configs <- list(
    "Informative" = c(0, 10, 3, 1),      # mu0, kappa0, alpha0, beta0
    "Weakly_Informative" = c(0, 1, 2, 2),
    "Diffuse" = c(0, 0.1, 1, 1)
  )

  results <- list()

  for (prior_name in names(prior_configs)) {
    cat(sprintf("\nTesting %s prior...\n", prior_name))

    prior_params <- prior_configs[[prior_name]]

    expr_list_prior <- list()

    # R implementation with specific prior
    expr_list_prior[[paste0("R_", prior_name)]] <- substitute({
      set_use_cpp(FALSE)
      mdobj <- MixingDistribution("normal", priorParameters = prior_params, "conjugate")
      dp <- DirichletProcessCreate(data, mdobj)
      dp <- Fit(dp, 50, progressBar = FALSE)
    }, list(prior_params = prior_params))

    # C++ implementation with specific prior
    if (exists("_dirichletprocess_run_mcmc_cpp")) {
      expr_list_prior[[paste0("Cpp_", prior_name)]] <- substitute({
        set_use_cpp(TRUE)
        mdobj <- MixingDistribution("normal", priorParameters = prior_params, "conjugate")
        dp <- DirichletProcessCreate(data, mdobj)
        dp <- Fit(dp, 50, progressBar = FALSE)
      }, list(prior_params = prior_params))
    }

    # Run benchmark for this prior configuration
    prior_result <- atime::atime(
      N = as.integer(10^seq(2, 3, by = 0.5)),
      setup = {
        data <- generate_test_data(N)
      },
      expr.list = expr_list_prior,
      seconds.limit = 5
    )

    results[[prior_name]] <- prior_result
    print(prior_result)
  }

  return(results)
}

# ==============================================================================
# Advanced atime Features
# ==============================================================================

#' Use atime_grid for comprehensive parameter sweep
run_grid_benchmark <- function() {

  cat("\n\n========================================\n")
  cat("Grid Benchmark (atime_grid)\n")
  cat("========================================\n\n")

  # Parameter grid
  param_grid <- atime::atime_grid(
    list(
      data_size = as.integer(10^seq(2, 3, by = 0.5)),
      n_iterations = c(50, 100, 200),
      n_clusters = c(2, 3, 5)
    ),
    expr = {
      data <- generate_test_data(data_size, k_clusters = n_clusters)

      # R implementation
      set_use_cpp(FALSE)
      dp_r <- DirichletProcessGaussian(data)
      dp_r <- Fit(dp_r, n_iterations, progressBar = FALSE)

      # C++ implementation (if available)
      if (exists("_dirichletprocess_run_mcmc_cpp")) {
        set_use_cpp(TRUE)
        dp_cpp <- DirichletProcessGaussian(data)
        dp_cpp <- Fit(dp_cpp, n_iterations, progressBar = FALSE)
      }
    }
  )

  print(param_grid)

  return(param_grid)
}

# ==============================================================================
# Main Execution
# ==============================================================================

if (interactive()) {
  cat("Normal Distribution atime Benchmark Suite\n")
  cat("=========================================\n\n")
  cat("Available benchmarks:\n")
  cat("  run_atime_benchmark()       - Main scaling comparison\n")
  cat("  benchmark_components()      - Component-level analysis\n")
  cat("  analyze_memory_scaling()    - Memory usage scaling\n")
  cat("  benchmark_prior_variations() - Different prior configurations\n")
  cat("  run_grid_benchmark()        - Comprehensive parameter sweep\n")
  cat("\nRun any function to start benchmarking!\n")

  # Quick test to verify setup
  cat("\nRunning quick verification test...\n")
  test_data <- generate_test_data(100)

  set_use_cpp(FALSE)
  dp_test <- DirichletProcessGaussian(test_data)
  dp_test <- Fit(dp_test, 10, progressBar = FALSE)
  cat(sprintf("✓ R implementation works (found %d clusters)\n", dp_test$numberClusters))

  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    set_use_cpp(TRUE)
    dp_test_cpp <- DirichletProcessGaussian(test_data)
    dp_test_cpp <- Fit(dp_test_cpp, 10, progressBar = FALSE)
    cat(sprintf("✓ C++ implementation works (found %d clusters)\n", dp_test_cpp$numberClusters))
  }

  cat("\nReady for benchmarking!\n")
}
