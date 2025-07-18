# benchmark-weibull-atime.R
# Benchmark Weibull Distribution using atime package
# Based on Neal (2000) and Escobar & West (1995) algorithms

library(dirichletprocess)
library(atime)
library(ggplot2)

# ==============================================================================
# Setup Functions
# ==============================================================================

#' Generate synthetic Weibull mixture data for benchmarking
#' Creates mixture of Weibull distributions with known clusters
generate_weibull_test_data <- function(n, k_clusters = 3, seed = 42) {
  set.seed(seed)

  # Equal sized clusters
  cluster_sizes <- rep(n %/% k_clusters, k_clusters)
  cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)

  # Well-separated Weibull parameters (shape, scale)
  # These create distinct shapes: decreasing hazard, constant hazard, increasing hazard
  shapes <- list(
    c(0.8, 2.0),   # Decreasing hazard rate
    c(2.0, 3.0),   # Increasing hazard rate
    c(1.0, 1.5)    # Constant hazard rate (exponential)
  )[1:k_clusters]

  data <- numeric(n)
  idx <- 1

  for (i in 1:k_clusters) {
    cluster_data <- rweibull(cluster_sizes[i],
                             shape = shapes[[i]][1],
                             scale = shapes[[i]][2])
    data[idx:(idx + cluster_sizes[i] - 1)] <- cluster_data
    idx <- idx + cluster_sizes[i]
  }

  return(data)
}

# ==============================================================================
# Main atime Benchmark
# ==============================================================================

#' Run atime benchmark comparing R and C++ implementations for Weibull
run_weibull_atime_benchmark <- function() {

  cat("========================================\n")
  cat("Weibull Distribution DP Benchmark (atime)\n")
  cat("========================================\n\n")

  # Check C++ availability for Weibull
  cpp_available <- exists("_dirichletprocess_run_mcmc_cpp") &&
    exists("can_use_cpp")

  if (!cpp_available) {
    cat("⚠️  C++ implementation not available for Weibull\n")
    cat("Only R implementation will be benchmarked\n\n")
  } else {
    # Test if C++ supports Weibull
    test_data <- rweibull(10, 2, 2)
    test_dp <- DirichletProcessWeibull(test_data, c(10, 2, 0.01), verbose = FALSE)
    cpp_supports_weibull <- tryCatch({
      can_use_cpp(test_dp)
    }, error = function(e) FALSE)

    if (!cpp_supports_weibull) {
      cat("⚠️  C++ implementation does not support Weibull distribution yet\n")
      cat("Only R implementation will be benchmarked\n\n")
      cpp_available <- FALSE
    } else {
      cat("✓ C++ implementation available for Weibull\n\n")
    }
  }

  # Define expression list for atime
  expr_list <- list()

  # R implementation expression
  expr_list$R_implementation <- quote({
    set_use_cpp(FALSE)
    set.seed(123)
    dp <- DirichletProcessWeibull(data, g0Priors, verbose = FALSE)
    dp <- Fit(dp, iterations, progressBar = FALSE)
  })

  # C++ implementation expression (if available)
  if (cpp_available) {
    expr_list$Cpp_implementation <- quote({
      set_use_cpp(TRUE)
      set.seed(123)
      dp <- DirichletProcessWeibull(data, g0Priors, verbose = FALSE)
      dp <- Fit(dp, iterations, progressBar = FALSE)
    })
  }

  # Run atime benchmark
  atime_result <- atime::atime(
    N = as.integer(10^seq(1.5, 3.5, by = 0.25)),  # 30 to 3000+ observations
    setup = {
      data <- generate_weibull_test_data(N, k_clusters = 3)
      iterations <- 100  # Fixed number of iterations
      g0Priors <- c(10, 2, 0.01)  # (phi, alpha0, beta0)
    },
    expr.list = expr_list,
    seconds.limit = 500,  # Stop if any expression takes more than 500 seconds
    verbose = TRUE
  )

  # Print summary
  print(atime_result)

  # Create and save plot
  p <- plot(atime_result)
  print(p)

  # Save results
  cat("\nSaving results...\n")
  save(atime_result, file = "atime_weibull_results.RData")
  ggsave("atime_weibull_benchmark.png", plot = p, width = 10, height = 8, dpi = 300)

  return(atime_result)
}

# ==============================================================================
# Component-level Benchmarking with atime
# ==============================================================================

#' Benchmark individual components for Weibull using atime
benchmark_weibull_components <- function() {

  cat("\n\n==========================================\n")
  cat("Weibull Component-level Benchmarking (atime)\n")
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
        mdObj <- WeibullMixtureCreate(c(10, 2, 0.01), mhStepSize = c(0.1, 0.1))
        # Create theta with proper array structure
        theta <- list(
          array(2.0, dim = c(1, 1, 1)),  # shape parameter
          array(3.0, dim = c(1, 1, 1))   # scale parameter
        )
        # Calculate likelihood for all data points
        for (i in 1:10) {
          lik <- Likelihood(mdObj, data, theta)
        }
      })
    } else if (comp == "PriorDraw") {
      expr_list_comp[["Prior_sampling"]] <- quote({
        mdObj <- WeibullMixtureCreate(c(10, 2, 0.01), mhStepSize = c(0.1, 0.1))
        # Draw multiple prior samples
        for (i in 1:10) {
          prior_samples <- PriorDraw(mdObj, 10)
        }
      })
    } else if (comp == "PosteriorDraw") {
      expr_list_comp[["Posterior_sampling"]] <- quote({
        mdObj <- WeibullMixtureCreate(c(10, 2, 0.01), mhStepSize = c(0.1, 0.1))
        # Draw posterior samples using MH
        post_samples <- PosteriorDraw(mdObj, matrix(data, ncol = 1), n = 10)
      })
    }

    # Run component benchmark
    comp_result <- atime::atime(
      N = as.integer(10^seq(2, 3.5, by = 0.5)),
      setup = {
        data <- generate_weibull_test_data(N)
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

#' Analyze memory scaling for Weibull distribution
analyze_weibull_memory_scaling <- function() {

  cat("\n\n========================================\n")
  cat("Weibull Memory Scaling Analysis (atime)\n")
  cat("========================================\n\n")

  # Focus on memory measurement
  expr_list_mem <- list()

  expr_list_mem$R_memory <- quote({
    set_use_cpp(FALSE)
    dp <- DirichletProcessWeibull(data, g0Priors, verbose = FALSE)
    dp <- Fit(dp, 50, progressBar = FALSE)
    # Force garbage collection to get accurate memory usage
    gc()
  })

  # Run memory-focused benchmark
  memory_result <- atime::atime(
    N = as.integer(10^seq(2, 4, by = 0.5)),  # Up to 10,000 observations
    setup = {
      data <- generate_weibull_test_data(N)
      g0Priors <- c(10, 2, 0.01)
    },
    expr.list = expr_list_mem,
    seconds.limit = 20
  )

  print(memory_result)
  plot(memory_result)

  return(memory_result)
}

# ==============================================================================
# Different Prior Configurations
# ==============================================================================

#' Benchmark with different prior configurations for Weibull
benchmark_weibull_prior_variations <- function() {

  cat("\n\n========================================\n")
  cat("Weibull Prior Configuration Benchmarks\n")
  cat("========================================\n\n")

  # Different prior configurations
  # c(phi, alpha0, beta0) - phi is upper bound for shape, alpha0/beta0 for scale prior
  prior_configs <- list(
    "Informative" = c(5, 3, 0.1),      # Tight bounds, informative
    "Weakly_Informative" = c(10, 2, 0.01),  # Default
    "Diffuse" = c(20, 1, 0.001)        # Wide bounds, less informative
  )

  results <- list()

  for (prior_name in names(prior_configs)) {
    cat(sprintf("\nTesting %s prior...\n", prior_name))

    prior_params <- prior_configs[[prior_name]]

    expr_list_prior <- list()

    # R implementation with specific prior
    expr_list_prior[[paste0("R_", prior_name)]] <- substitute({
      set_use_cpp(FALSE)
      set.seed(123)
      dp <- DirichletProcessWeibull(data, prior_params, verbose = FALSE)
      dp <- Fit(dp, 50, progressBar = FALSE)
    }, list(prior_params = prior_params))

    # Run benchmark for this prior configuration
    prior_result <- atime::atime(
      N = as.integer(10^seq(2, 3, by = 0.5)),
      setup = {
        data <- generate_weibull_test_data(N)
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
# Quick Performance Test
# ==============================================================================

#' Quick performance test at fixed sizes
test_weibull_performance <- function() {

  cat("\n==========================================\n")
  cat("Quick Weibull Performance Test\n")
  cat("==========================================\n\n")

  n_values <- c(50, 100, 200, 500)
  results <- data.frame()

  for (n in n_values) {
    cat(sprintf("Testing N=%d...\n", n))

    data <- generate_weibull_test_data(n, k_clusters = 3)
    g0Priors <- c(10, 2, 0.01)

    # Time R implementation
    set_use_cpp(FALSE)
    time_r <- system.time({
      dp <- DirichletProcessWeibull(data, g0Priors, verbose = FALSE)
      dp <- Fit(dp, 50, progressBar = FALSE)
    })[3]

    # Count clusters found
    n_clusters_r <- dp$numberClusters

    results <- rbind(results, data.frame(
      N = n,
      R_time = time_r,
      Clusters_found = n_clusters_r
    ))
  }

  print(results)

  return(results)
}

# ==============================================================================
# Comparison with Different Numbers of Clusters
# ==============================================================================

#' Benchmark performance vs number of true clusters
benchmark_weibull_cluster_scaling <- function() {

  cat("\n\n========================================\n")
  cat("Weibull Cluster Scaling Analysis\n")
  cat("========================================\n\n")

  cluster_counts <- c(2, 3, 5)
  results <- list()

  for (k in cluster_counts) {
    cat(sprintf("\nBenchmarking with %d clusters...\n", k))

    expr_list_k <- list()

    expr_list_k[[paste0("R_", k, "_clusters")]] <- quote({
      set_use_cpp(FALSE)
      set.seed(123)
      dp <- DirichletProcessWeibull(data, g0Priors, verbose = FALSE)
      dp <- Fit(dp, 100, progressBar = FALSE)
    })

    # Run benchmark
    k_result <- atime::atime(
      N = as.integer(c(100, 200, 500, 1000)),
      setup = {
        data <- generate_weibull_test_data(N, k_clusters = k)
        g0Priors <- c(10, 2, 0.01)
      },
      expr.list = expr_list_k,
      seconds.limit = 10
    )

    results[[paste0(k, "_clusters")]] <- k_result
    print(k_result)
  }

  return(results)
}

# ==============================================================================
# Main Execution
# ==============================================================================

if (interactive()) {
  cat("Weibull Distribution atime Benchmark Suite\n")
  cat("==========================================\n\n")
  cat("Available benchmarks:\n")
  cat("  run_weibull_atime_benchmark()      - Main scaling comparison\n")
  cat("  benchmark_weibull_components()     - Component-level analysis\n")
  cat("  analyze_weibull_memory_scaling()   - Memory usage scaling\n")
  cat("  benchmark_weibull_prior_variations() - Different prior configurations\n")
  cat("  test_weibull_performance()         - Quick performance test\n")
  cat("  benchmark_weibull_cluster_scaling() - Performance vs number of clusters\n")
  cat("\nRun any function to start benchmarking!\n")

  # Quick test to verify setup
  cat("\nRunning quick verification test...\n")
  test_data <- generate_weibull_test_data(100)

  set_use_cpp(FALSE)
  dp_test <- DirichletProcessWeibull(test_data, c(10, 2, 0.01), verbose = FALSE)
  dp_test <- Fit(dp_test, 10, progressBar = FALSE)
  cat(sprintf("✓ R implementation works (found %d clusters)\n", dp_test$numberClusters))

  cat("\nReady for benchmarking!\n")
}
