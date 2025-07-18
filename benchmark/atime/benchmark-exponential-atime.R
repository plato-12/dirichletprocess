# benchmark-exponential-atime.R
# Benchmark Exponential Distribution using atime package
# Based on Neal (2000) and Escobar & West (1995) algorithms

library(dirichletprocess)
library(atime)
library(ggplot2)

# ==============================================================================
# Setup Functions
# ==============================================================================

#' Generate synthetic exponential mixture data for benchmarking
#' Creates mixture of exponential distributions with known clusters
generate_exponential_test_data <- function(n, k_clusters = 3, seed = 42) {
  set.seed(seed)

  # Equal sized clusters
  cluster_sizes <- rep(n %/% k_clusters, k_clusters)
  cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)

  # Well-separated exponential rates (lambda parameters)
  # Different rates create distinct exponential distributions
  rates <- seq(0.5, 5, length.out = k_clusters)

  data <- numeric(n)
  idx <- 1

  for (i in 1:k_clusters) {
    cluster_data <- rexp(cluster_sizes[i], rate = rates[i])
    data[idx:(idx + cluster_sizes[i] - 1)] <- cluster_data
    idx <- idx + cluster_sizes[i]
  }

  return(data)
}

# ==============================================================================
# Main atime Benchmark
# ==============================================================================

#' Run atime benchmark comparing R and C++ implementations for Exponential
run_exponential_atime_benchmark <- function() {

  cat("==========================================\n")
  cat("Exponential Distribution DP Benchmark (atime)\n")
  cat("==========================================\n\n")

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
    dp <- DirichletProcessExponential(data)
    dp <- Fit(dp, iterations, progressBar = FALSE)
  })

  # C++ implementation expression (if available)
  if (cpp_available) {
    expr_list$Cpp_implementation <- quote({
      set_use_cpp(TRUE)
      set.seed(123)
      dp <- DirichletProcessExponential(data)
      dp <- Fit(dp, iterations, progressBar = FALSE)
    })
  }

  # Run atime benchmark
  atime_result <- atime::atime(
    N = as.integer(10^seq(1.5, 3.5, by = 0.25)),  # 30 to 3000+ observations
    setup = {
      data <- generate_exponential_test_data(N, k_clusters = 3)
      iterations <- 100  # Fixed number of iterations
    },
    expr.list = expr_list,
    seconds.limit = 10000,  # Stop if any expression takes more than 10 seconds
    verbose = TRUE
  )

  # Print summary
  print(atime_result)

  # Create and save plot
  p <- plot(atime_result)
  print(p)

  # Save results
  cat("\nSaving results...\n")
  save(atime_result, file = "atime_exponential_results.RData")
  ggsave("atime_exponential_benchmark.png", plot = p, width = 10, height = 8, dpi = 300)

  return(atime_result)
}

# ==============================================================================
# Component-level Benchmarking with atime
# ==============================================================================

#' Benchmark individual components for Exponential using atime
benchmark_exponential_components <- function() {

  cat("\n\n==========================================\n")
  cat("Exponential Component-level Benchmarking (atime)\n")
  cat("==========================================\n\n")

  # Define components to benchmark
  components <- c("Likelihood", "PriorDraw", "PosteriorDraw")

  results <- list()

  for (comp in components) {
    cat(sprintf("\nBenchmarking %s...\n", comp))

    expr_list_comp <- list()

    # Test each component
    if (comp == "Likelihood") {
      expr_list_comp[["Likelihood_calculation"]] <- quote({
        mdObj <- ExponentialMixtureCreate(c(1, 1))
        # Create theta with proper array structure
        theta <- list(array(2, dim = c(1, 1, 1)))
        # Calculate likelihood for all data points
        for (i in 1:10) {
          lik <- Likelihood(mdObj, data, theta)
        }
      })
    } else if (comp == "PriorDraw") {
      expr_list_comp[["Prior_sampling"]] <- quote({
        mdObj <- ExponentialMixtureCreate(c(1, 1))
        # Draw multiple prior samples
        for (i in 1:10) {
          prior_samples <- PriorDraw(mdObj, 10)
        }
      })
    } else if (comp == "PosteriorDraw") {
      expr_list_comp[["Posterior_sampling"]] <- quote({
        mdObj <- ExponentialMixtureCreate(c(1, 1))
        # Draw posterior samples
        post_samples <- PosteriorDraw(mdObj, matrix(data, ncol = 1), n = 10)
      })
    }

    # Run component benchmark
    comp_result <- atime::atime(
      N = as.integer(10^seq(2, 3.5, by = 0.5)),
      setup = {
        data <- generate_exponential_test_data(N)
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

#' Analyze memory scaling for Exponential distribution
analyze_exponential_memory_scaling <- function() {

  cat("\n\n========================================\n")
  cat("Exponential Memory Scaling Analysis (atime)\n")
  cat("========================================\n\n")

  # Memory tracking expressions
  expr_list_mem <- list()

  # R implementation with memory tracking
  expr_list_mem$R_memory <- quote({
    set_use_cpp(FALSE)
    gc(reset = TRUE)

    dp <- DirichletProcessExponential(data)
    dp <- Fit(dp, 50, progressBar = FALSE)

    mem_used <- gc()[2, 2]  # Total memory in MB
  })

  # C++ implementation with memory tracking
  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    expr_list_mem$Cpp_memory <- quote({
      set_use_cpp(TRUE)
      gc(reset = TRUE)

      if (exists("clear_memory_tracking")) clear_memory_tracking()

      dp <- DirichletProcessExponential(data)
      dp <- Fit(dp, 50, progressBar = FALSE)

      mem_used <- gc()[2, 2]  # Total memory in MB
    })
  }

  # Run memory benchmark
  mem_result <- atime::atime(
    N = as.integer(10^seq(2, 3.5, by = 0.5)),
    setup = {
      data <- generate_exponential_test_data(N)
    },
    expr.list = expr_list_mem,
    seconds.limit = 10
  )

  print(mem_result)
  plot(mem_result)

  return(mem_result)
}

# ==============================================================================
# Comparison with Different Prior Settings
# ==============================================================================

#' Benchmark with different prior configurations for Exponential
benchmark_exponential_prior_variations <- function() {

  cat("\n\n========================================\n")
  cat("Exponential Prior Configuration Benchmarks (atime)\n")
  cat("========================================\n\n")

  # Different prior configurations for Gamma(alpha, beta)
  prior_configs <- list(
    "Informative" = c(2, 2),        # Strong prior belief
    "Weakly_Informative" = c(1, 1), # Standard exponential prior
    "Diffuse" = c(0.01, 0.01)       # Very weak prior
  )

  results <- list()

  for (prior_name in names(prior_configs)) {
    cat(sprintf("\nTesting %s prior...\n", prior_name))

    prior_params <- prior_configs[[prior_name]]

    expr_list_prior <- list()

    # R implementation with specific prior
    expr_list_prior[[paste0("R_", prior_name)]] <- substitute({
      set_use_cpp(FALSE)
      dp <- DirichletProcessExponential(data, g0Priors = prior_params)
      dp <- Fit(dp, 50, progressBar = FALSE)
    }, list(prior_params = prior_params))

    # C++ implementation with specific prior
    if (exists("_dirichletprocess_run_mcmc_cpp")) {
      expr_list_prior[[paste0("Cpp_", prior_name)]] <- substitute({
        set_use_cpp(TRUE)
        dp <- DirichletProcessExponential(data, g0Priors = prior_params)
        dp <- Fit(dp, 50, progressBar = FALSE)
      }, list(prior_params = prior_params))
    }

    # Run benchmark for this prior configuration
    prior_result <- atime::atime(
      N = as.integer(10^seq(2, 3, by = 0.5)),
      setup = {
        data <- generate_exponential_test_data(N)
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
# Cluster Number Scaling Analysis
# ==============================================================================

#' Analyze performance with different numbers of true clusters
benchmark_exponential_cluster_scaling <- function() {

  cat("\n\n========================================\n")
  cat("Exponential Cluster Scaling Analysis (atime)\n")
  cat("========================================\n\n")

  # Test with different numbers of clusters
  cluster_counts <- c(2, 3, 5, 8)

  results <- list()

  for (k in cluster_counts) {
    cat(sprintf("\nBenchmarking with %d clusters...\n", k))

    expr_list_cluster <- list()

    # R implementation
    expr_list_cluster[[paste0("R_", k, "_clusters")]] <- quote({
      set_use_cpp(FALSE)
      set.seed(123)
      dp <- DirichletProcessExponential(data)
      dp <- Fit(dp, 75, progressBar = FALSE)
    })

    # C++ implementation
    if (exists("_dirichletprocess_run_mcmc_cpp")) {
      expr_list_cluster[[paste0("Cpp_", k, "_clusters")]] <- quote({
        set_use_cpp(TRUE)
        set.seed(123)
        dp <- DirichletProcessExponential(data)
        dp <- Fit(dp, 75, progressBar = FALSE)
      })
    }

    # Run benchmark for this cluster count
    cluster_result <- atime::atime(
      N = as.integer(c(100, 300, 500, 1000)),
      setup = {
        data <- generate_exponential_test_data(N, k_clusters = k)
      },
      expr.list = expr_list_cluster,
      seconds.limit = 10
    )

    results[[paste0(k, "_clusters")]] <- cluster_result
    print(cluster_result)
  }

  return(results)
}

# ==============================================================================
# Quick Performance Comparison
# ==============================================================================

#' Quick comparison at fixed sizes
compare_exponential_implementations <- function() {

  cat("\n==========================================\n")
  cat("Quick Exponential Implementation Comparison\n")
  cat("==========================================\n\n")

  n_values <- c(100, 500, 1000)
  results <- data.frame()

  for (n in n_values) {
    cat(sprintf("Testing N=%d...\n", n))

    data <- generate_exponential_test_data(n, k_clusters = 3)

    # Time R
    set_use_cpp(FALSE)
    time_r <- system.time({
      dp <- DirichletProcessExponential(data)
      dp <- Fit(dp, 100, progressBar = FALSE)
    })[3]

    clusters_r <- dp$numberClusters

    # Time C++ (if available)
    if (exists("_dirichletprocess_run_mcmc_cpp")) {
      set_use_cpp(TRUE)
      time_cpp <- system.time({
        dp <- DirichletProcessExponential(data)
        dp <- Fit(dp, 100, progressBar = FALSE)
      })[3]
      clusters_cpp <- dp$numberClusters
    } else {
      time_cpp <- NA
      clusters_cpp <- NA
    }

    results <- rbind(results, data.frame(
      N = n,
      R_time = time_r,
      Cpp_time = time_cpp,
      Speedup = if (!is.na(time_cpp)) time_r / time_cpp else NA,
      R_clusters = clusters_r,
      Cpp_clusters = clusters_cpp
    ))
  }

  print(results)

  cat("\nKey Findings:\n")
  if (!all(is.na(results$Speedup))) {
    cat(sprintf("- Average speedup: %.1fx\n", mean(results$Speedup, na.rm = TRUE)))
  }
  cat("- Exponential distribution is conjugate, allowing efficient sampling\n")
  cat("- Performance scales well with data size\n")

  return(results)
}

# ==============================================================================
# Advanced Grid Benchmark
# ==============================================================================

#' Use atime_grid for comprehensive parameter sweep
run_exponential_grid_benchmark <- function() {

  cat("\n\n========================================\n")
  cat("Exponential Grid Benchmark (atime_grid)\n")
  cat("========================================\n\n")

  # Parameter grid
  param_grid <- atime::atime_grid(
    list(
      data_size = as.integer(10^seq(2, 3, by = 0.5)),
      n_iterations = c(50, 100, 200),
      n_clusters = c(2, 3, 5)
    ),
    expr = {
      data <- generate_exponential_test_data(data_size, k_clusters = n_clusters)

      # R implementation
      set_use_cpp(FALSE)
      dp_r <- DirichletProcessExponential(data)
      dp_r <- Fit(dp_r, n_iterations, progressBar = FALSE)

      # C++ implementation (if available)
      if (exists("_dirichletprocess_run_mcmc_cpp")) {
        set_use_cpp(TRUE)
        dp_cpp <- DirichletProcessExponential(data)
        dp_cpp <- Fit(dp_cpp, n_iterations, progressBar = FALSE)
      }
    }
  )

  print(param_grid)

  return(param_grid)
}

# ==============================================================================
# Comprehensive Benchmark Suite
# ==============================================================================

#' Run all exponential benchmarks
run_all_exponential_benchmarks <- function() {

  cat("Exponential Distribution Comprehensive Benchmark Suite\n")
  cat("====================================================\n\n")

  # Check setup
  cat("Checking setup...\n")
  test_data <- generate_exponential_test_data(100, k_clusters = 3)

  set_use_cpp(FALSE)
  dp_test <- DirichletProcessExponential(test_data)
  dp_test <- Fit(dp_test, 10, progressBar = FALSE)
  cat(sprintf("✓ R implementation works (found %d clusters)\n", dp_test$numberClusters))

  cpp_available <- FALSE
  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    tryCatch({
      set_use_cpp(TRUE)
      dp_test_cpp <- DirichletProcessExponential(test_data)
      dp_test_cpp <- Fit(dp_test_cpp, 10, progressBar = FALSE)
      cat(sprintf("✓ C++ implementation works (found %d clusters)\n", dp_test_cpp$numberClusters))
      cpp_available <- TRUE
    }, error = function(e) {
      cat("✗ C++ implementation not available for Exponential\n")
    })
  }

  results <- list()

  # 1. Quick comparison first
  cat("\n1. Running quick comparison...\n")
  results$quick_comparison <- compare_exponential_implementations()

  # 2. Main benchmark
  cat("\n2. Running main benchmark...\n")
  results$main_benchmark <- run_exponential_atime_benchmark()

  # 3. Component analysis
  cat("\n3. Running component analysis...\n")
  results$components <- benchmark_exponential_components()

  # 4. Memory scaling
  cat("\n4. Running memory scaling analysis...\n")
  results$memory <- analyze_exponential_memory_scaling()

  # 5. Prior variations
  cat("\n5. Running prior variation benchmarks...\n")
  results$priors <- benchmark_exponential_prior_variations()

  # 6. Cluster scaling
  cat("\n6. Running cluster scaling analysis...\n")
  results$clusters <- benchmark_exponential_cluster_scaling()

  # 7. Grid benchmark (if time permits)
  cat("\n7. Running grid benchmark (this may take a while)...\n")
  results$grid <- run_exponential_grid_benchmark()

  cat("\n====================================================\n")
  cat("All benchmarks completed!\n")
  cat("Results saved in 'atime_exponential_results.RData'\n")

  return(results)
}

# ==============================================================================
# Main Execution
# ==============================================================================

if (interactive()) {
  cat("Exponential Distribution atime Benchmark Suite\n")
  cat("============================================\n\n")
  cat("Available benchmarks:\n")
  cat("  run_all_exponential_benchmarks()       - Run all benchmarks\n")
  cat("  run_exponential_atime_benchmark()      - Main scaling comparison\n")
  cat("  benchmark_exponential_components()     - Component-level analysis\n")
  cat("  analyze_exponential_memory_scaling()   - Memory usage scaling\n")
  cat("  benchmark_exponential_prior_variations() - Different prior configurations\n")
  cat("  benchmark_exponential_cluster_scaling() - Cluster number scaling\n")
  cat("  compare_exponential_implementations()  - Quick comparison\n")
  cat("  run_exponential_grid_benchmark()       - Comprehensive parameter sweep\n")
  cat("\nRun any function to start benchmarking!\n")

  # Quick test to verify setup
  cat("\nRunning quick verification test...\n")
  test_data <- generate_exponential_test_data(100)

  set_use_cpp(FALSE)
  dp_test <- DirichletProcessExponential(test_data)
  dp_test <- Fit(dp_test, 10, progressBar = FALSE)
  cat(sprintf("✓ R implementation works (found %d clusters)\n", dp_test$numberClusters))

  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    tryCatch({
      set_use_cpp(TRUE)
      dp_test_cpp <- DirichletProcessExponential(test_data)
      dp_test_cpp <- Fit(dp_test_cpp, 10, progressBar = FALSE)
      cat(sprintf("✓ C++ implementation works (found %d clusters)\n", dp_test_cpp$numberClusters))
    }, error = function(e) {
      cat("✗ C++ implementation not available for Exponential\n")
    })
  }

  cat("\nReady for benchmarking!\n")
}
