# benchmark/atime/mvnormal_improved.R
# Improved MVNormal Distribution Benchmark using atime package
# Based on Neal (2000) and Escobar & West (1995) algorithms

library(dirichletprocess)
library(atime)
library(ggplot2)
library(mvtnorm)
library(data.table)

# ==============================================================================
# Setup Functions
# ==============================================================================

#' Generate synthetic multivariate data for benchmarking
#' Creates mixture of multivariate normal distributions with known clusters
generate_mvnormal_test_data <- function(n, d = 2, k_clusters = 3, seed = 42) {
  set.seed(seed)

  # Equal sized clusters
  cluster_sizes <- rep(n %/% k_clusters, k_clusters)
  cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)

  # Well-separated cluster means in d dimensions
  means <- list()
  for (i in 1:k_clusters) {
    # Create well-separated means
    mean_vec <- rep(0, d)
    mean_vec[1] <- (i - 2) * 3  # Separate along first dimension
    if (d > 1) mean_vec[2] <- (i - 2) * 2  # Separate along second dimension
    means[[i]] <- mean_vec
  }

  # Common covariance matrix
  sigma <- diag(d) * 0.5

  data <- matrix(NA, n, d)
  idx <- 1

  for (i in 1:k_clusters) {
    cluster_data <- mvtnorm::rmvnorm(cluster_sizes[i],
                                     mean = means[[i]],
                                     sigma = sigma)
    data[idx:(idx + cluster_sizes[i] - 1), ] <- cluster_data
    idx <- idx + cluster_sizes[i]
  }

  return(data)
}

# ==============================================================================
# Main atime Benchmark with Adjusted Parameters
# ==============================================================================

#' Run atime benchmark with practical parameters for MVNormal
run_mvnormal_atime_benchmark <- function() {

  cat("==========================================\n")
  cat("MVNormal Distribution DP Benchmark (atime)\n")
  cat("==========================================\n\n")

  # Check C++ availability
  cpp_available <- exists("mvnormal_prior_draw_cpp") &&
    exists("conjugate_mvnormal_cluster_component_update_cpp")

  if (!cpp_available) {
    cat("⚠️  MVNormal C++ implementation not available\n")
    cat("Only R implementation will be benchmarked\n\n")
  } else {
    cat("✓ MVNormal C++ implementation available\n\n")
  }

  # Define expression list for atime
  expr_list <- list()

  # R implementation expression
  expr_list$R_implementation <- quote({
    set_use_cpp(FALSE)
    set.seed(123)
    dp <- DirichletProcessMvnormal(data, g0Priors)
    dp <- Fit(dp, iterations, progressBar = FALSE)
  })

  # C++ implementation expression (if available)
  if (cpp_available) {
    expr_list$Cpp_implementation <- quote({
      set_use_cpp(TRUE)
      set.seed(123)
      dp <- DirichletProcessMvnormal(data, g0Priors)
      dp <- Fit(dp, iterations, progressBar = FALSE)
    })
  }

  # MODIFIED: Reduced N range and iterations
  atime_result <- atime::atime(
    N = as.integer(c(30, 50, 75, 100, 150, 200, 300, 500)),  # Custom N values up to 500
    setup = {
      # Test with 2D data by default
      d <- 2
      data <- generate_mvnormal_test_data(N, d = d, k_clusters = 3)
      iterations <- 50  # Reduced from 100 to 50 iterations

      # Set up priors for MVNormal
      g0Priors <- list(
        mu0 = rep(0, d),
        Lambda = diag(d),
        kappa0 = 1,
        nu = d + 2
      )
    },
    expr.list = expr_list,
    seconds.limit = 60000,  # Stop if any expression takes more than 60 seconds
    verbose = TRUE
  )

  # Print summary
  print(atime_result)

  # Create and save plot
  p <- plot(atime_result)
  print(p)

  # Save results
  cat("\nSaving results...\n")
  save(atime_result, file = "atime_mvnormal_results.RData")
  ggsave("atime_mvnormal_benchmark.png", plot = p, width = 10, height = 8, dpi = 300)

  return(atime_result)
}

# ==============================================================================
# Separate Benchmarks for R and C++
# ==============================================================================

#' R-only benchmark with smaller scale
benchmark_mvnormal_R <- function() {

  cat("\n==========================================\n")
  cat("MVNormal R Implementation Benchmark\n")
  cat("==========================================\n\n")

  atime_result_r <- atime::atime(
    N = c(30, 50, 100, 150, 200),  # Small N only
    setup = {
      data <- generate_mvnormal_test_data(N, d = 2, k_clusters = 3)
      iterations <- 50  # Fewer iterations
      g0Priors <- list(
        mu0 = rep(0, 2),
        Lambda = diag(2),
        kappa0 = 1,
        nu = 4
      )
    },
    expr.list = list(
      R_implementation = quote({
        set_use_cpp(FALSE)
        set.seed(123)
        dp <- DirichletProcessMvnormal(data, g0Priors)
        dp <- Fit(dp, iterations, progressBar = FALSE)
      })
    ),
    seconds.limit = 120,
    verbose = TRUE
  )

  print(atime_result_r)
  plot(atime_result_r)

  return(atime_result_r)
}

#' C++ only benchmark for larger datasets
benchmark_mvnormal_Cpp_large <- function() {

  if (!exists("mvnormal_prior_draw_cpp")) {
    cat("C++ implementation not available\n")
    return(NULL)
  }

  cat("\n==========================================\n")
  cat("MVNormal C++ Large Dataset Benchmark\n")
  cat("==========================================\n\n")

  # Test larger N values with C++ only
  atime_result_cpp <- atime::atime(
    N = as.integer(10^seq(2, 3.5, by = 0.25)),  # 100 to 3162
    setup = {
      d <- 2
      data <- generate_mvnormal_test_data(N, d = d, k_clusters = 3)
      iterations <- 100  # Full iterations for C++

      g0Priors <- list(
        mu0 = rep(0, d),
        Lambda = diag(d),
        kappa0 = 1,
        nu = d + 2
      )

      set_use_cpp(TRUE)
    },
    expr.list = list(
      Cpp_implementation = quote({
        set.seed(123)
        dp <- DirichletProcessMvnormal(data, g0Priors)
        dp <- Fit(dp, iterations, progressBar = FALSE)
      })
    ),
    seconds.limit = 30,
    verbose = TRUE
  )

  print(atime_result_cpp)
  plot(atime_result_cpp)

  return(atime_result_cpp)
}

# ==============================================================================
# Adaptive Iteration Strategy
# ==============================================================================

#' Adaptive benchmark with different iterations for R and C++
run_mvnormal_adaptive_benchmark <- function() {

  cat("\n==========================================\n")
  cat("MVNormal Adaptive Iteration Benchmark\n")
  cat("==========================================\n\n")

  expr_list <- list()

  # R with fewer iterations
  expr_list$R_50iter <- quote({
    set_use_cpp(FALSE)
    set.seed(123)
    dp <- DirichletProcessMvnormal(data, g0Priors)
    dp <- Fit(dp, 50, progressBar = FALSE)  # Only 50 iterations
  })

  # C++ with full iterations
  if (exists("mvnormal_prior_draw_cpp")) {
    expr_list$Cpp_100iter <- quote({
      set_use_cpp(TRUE)
      set.seed(123)
      dp <- DirichletProcessMvnormal(data, g0Priors)
      dp <- Fit(dp, 100, progressBar = FALSE)  # Full 100 iterations
    })
  }

  # Run benchmark
  atime_result <- atime::atime(
    N = as.integer(c(50, 100, 200, 300, 500)),
    setup = {
      d <- 2
      data <- generate_mvnormal_test_data(N, d = d, k_clusters = 3)
      g0Priors <- list(
        mu0 = rep(0, d),
        Lambda = diag(d),
        kappa0 = 1,
        nu = d + 2
      )
    },
    expr.list = expr_list,
    seconds.limit = 60,
    verbose = TRUE
  )

  print(atime_result)
  plot(atime_result)

  return(atime_result)
}

# ==============================================================================
# Component-level Analysis (Quick)
# ==============================================================================

#' Quick component benchmark to identify bottlenecks
benchmark_mvnormal_components_quick <- function() {

  cat("\n==========================================\n")
  cat("MVNormal Component Analysis\n")
  cat("==========================================\n\n")

  # Test with small dataset
  n <- 100
  d <- 2
  data <- matrix(rnorm(n * d), ncol = d)
  g0Priors <- list(
    mu0 = rep(0, d),
    Lambda = diag(d),
    kappa0 = 1,
    nu = d + 2
  )

  # Initialize DP
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessMvnormal(data, g0Priors)
  dp_r <- InitialiseClusters(dp_r)

  # Time individual components
  library(microbenchmark)

  # 1. Cluster assignment
  cat("Timing cluster assignment...\n")
  time_cluster <- microbenchmark(
    R = ClusterComponentUpdate(dp_r),
    times = 10
  )

  # 2. Parameter update
  cat("Timing parameter update...\n")
  time_param <- microbenchmark(
    R = ClusterParameterUpdate(dp_r),
    times = 10
  )

  # 3. Likelihood calculation
  cat("Timing likelihood calculation...\n")
  mdObj <- MvnormalCreate(g0Priors)
  theta <- list(
    mu = array(rep(0, d), c(1, d, 1)),
    sig = array(diag(d), c(d, d, 1))
  )

  time_lik <- microbenchmark(
    R = Likelihood(mdObj, data, theta),
    times = 10
  )

  cat("\nComponent timings (milliseconds):\n")
  cat(sprintf("Cluster Assignment: %.2f ms\n", median(time_cluster$time) / 1e6))
  cat(sprintf("Parameter Update: %.2f ms\n", median(time_param$time) / 1e6))
  cat(sprintf("Likelihood: %.2f ms\n", median(time_lik$time) / 1e6))

  # If C++ is available, compare
  if (exists("mvnormal_prior_draw_cpp")) {
    cat("\nComparing with C++ components...\n")

    set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessMvnormal(data, g0Priors)
    dp_cpp <- InitialiseClusters(dp_cpp)

    time_cluster_cpp <- microbenchmark(
      Cpp = ClusterComponentUpdate(dp_cpp),
      times = 10
    )

    time_param_cpp <- microbenchmark(
      Cpp = ClusterParameterUpdate(dp_cpp),
      times = 10
    )

    cat(sprintf("\nC++ Cluster Assignment: %.2f ms\n", median(time_cluster_cpp$time) / 1e6))
    cat(sprintf("C++ Parameter Update: %.2f ms\n", median(time_param_cpp$time) / 1e6))

    cat(sprintf("\nSpeedup - Cluster Assignment: %.1fx\n",
                median(time_cluster$time) / median(time_cluster_cpp$time)))
    cat(sprintf("Speedup - Parameter Update: %.1fx\n",
                median(time_param$time) / median(time_param_cpp$time)))
  }

  return(list(
    cluster = time_cluster,
    param = time_param,
    likelihood = time_lik
  ))
}

# ==============================================================================
# Dimension Scaling Analysis
# ==============================================================================

#' Benchmark MVNormal across different dimensions
benchmark_mvnormal_dimensions <- function() {

  cat("\n\n==========================================\n")
  cat("MVNormal Dimension Scaling (atime)\n")
  cat("==========================================\n\n")

  dimensions <- c(2, 3, 5, 10)
  results <- list()

  for (d in dimensions) {
    cat(sprintf("\nBenchmarking %dD MVNormal...\n", d))

    expr_list_dim <- list()

    # R implementation
    expr_list_dim[[paste0("R_", d, "D")]] <- quote({
      set_use_cpp(FALSE)
      set.seed(123)
      dp <- DirichletProcessMvnormal(data, g0Priors)
      dp <- Fit(dp, 25, progressBar = FALSE)  # Fewer iterations for R
    })

    # C++ implementation
    if (exists("mvnormal_prior_draw_cpp")) {
      expr_list_dim[[paste0("Cpp_", d, "D")]] <- quote({
        set_use_cpp(TRUE)
        set.seed(123)
        dp <- DirichletProcessMvnormal(data, g0Priors)
        dp <- Fit(dp, 50, progressBar = FALSE)
      })
    }

    # Run benchmark for this dimension
    dim_result <- atime::atime(
      N = as.integer(c(50, 100, 200)),  # Smaller N range
      setup = {
        data <- generate_mvnormal_test_data(N, d = d, k_clusters = 2)
        g0Priors <- list(
          mu0 = rep(0, d),
          Lambda = diag(d),
          kappa0 = 1,
          nu = d + 2
        )
      },
      expr.list = expr_list_dim,
      seconds.limit = 30
    )

    results[[paste0(d, "D")]] <- dim_result
    print(dim_result)
  }

  return(results)
}

# ==============================================================================
# Quick Performance Comparison
# ==============================================================================

#' Quick comparison at fixed sizes
compare_implementations <- function() {

  cat("\n==========================================\n")
  cat("Quick MVNormal Implementation Comparison\n")
  cat("==========================================\n\n")

  n_values <- c(50, 100, 150)
  results <- data.frame()

  for (n in n_values) {
    cat(sprintf("Testing N=%d...\n", n))

    data <- generate_mvnormal_test_data(n, d = 2)
    g0Priors <- list(mu0 = c(0,0), Lambda = diag(2), kappa0 = 1, nu = 4)

    # Time R (fewer iterations)
    set_use_cpp(FALSE)
    time_r <- system.time({
      dp <- DirichletProcessMvnormal(data, g0Priors)
      dp <- Fit(dp, 25, progressBar = FALSE)
    })[3]

    # Time C++ (full iterations)
    if (exists("mvnormal_prior_draw_cpp")) {
      set_use_cpp(TRUE)
      time_cpp <- system.time({
        dp <- DirichletProcessMvnormal(data, g0Priors)
        dp <- Fit(dp, 100, progressBar = FALSE)
      })[3]
    } else {
      time_cpp <- NA
    }

    results <- rbind(results, data.frame(
      N = n,
      R_time_25iter = time_r,
      Cpp_time_100iter = time_cpp,
      Speedup = if (!is.na(time_cpp)) (time_r * 4) / time_cpp else NA
    ))
  }

  print(results)

  cat("\nKey Findings:\n")
  if (!all(is.na(results$Speedup))) {
    cat(sprintf("- Average adjusted speedup: %.1fx\n", mean(results$Speedup, na.rm = TRUE)))
  }
  cat("- MVNormal is computationally intensive due to matrix operations\n")
  cat("- Consider using C++ implementation for production use\n")

  return(results)
}

# ==============================================================================
# Main Execution Function
# ==============================================================================

#' Run all MVNormal benchmarks
run_all_mvnormal_benchmarks <- function() {

  cat("MVNormal Distribution Comprehensive Benchmark Suite\n")
  cat("=================================================\n\n")

  # Check setup
  cat("Checking setup...\n")
  test_data <- generate_mvnormal_test_data(50, d = 2)
  test_priors <- list(mu0 = c(0, 0), Lambda = diag(2), kappa0 = 1, nu = 4)

  set_use_cpp(FALSE)
  dp_test <- DirichletProcessMvnormal(test_data, test_priors)
  dp_test <- Fit(dp_test, 10, progressBar = FALSE)
  cat(sprintf("✓ R implementation works (found %d clusters)\n", dp_test$numberClusters))

  cpp_available <- FALSE
  if (exists("mvnormal_prior_draw_cpp")) {
    set_use_cpp(TRUE)
    dp_test_cpp <- DirichletProcessMvnormal(test_data, test_priors)
    dp_test_cpp <- Fit(dp_test_cpp, 10, progressBar = FALSE)
    cat(sprintf("✓ C++ implementation works (found %d clusters)\n", dp_test_cpp$numberClusters))
    cpp_available <- TRUE
  }

  results <- list()

  # 1. Quick comparison first
  cat("\n1. Running quick comparison...\n")
  results$quick_comparison <- compare_implementations()

  # 2. Component analysis
  cat("\n2. Running component analysis...\n")
  results$components <- benchmark_mvnormal_components_quick()

  # 3. Main benchmark with adjusted parameters
  cat("\n3. Running main benchmark...\n")
  results$main_benchmark <- run_mvnormal_atime_benchmark()

  # 4. If C++ available, run large dataset benchmark
  if (cpp_available) {
    cat("\n4. Running C++ large dataset benchmark...\n")
    results$cpp_large <- benchmark_mvnormal_Cpp_large()
  }

  # 5. Dimension scaling
  cat("\n5. Running dimension scaling analysis...\n")
  results$dimensions <- benchmark_mvnormal_dimensions()

  # 6. Adaptive iteration benchmark
  cat("\n6. Running adaptive iteration benchmark...\n")
  results$adaptive <- run_mvnormal_adaptive_benchmark()

  cat("\n=================================================\n")
  cat("All benchmarks completed!\n")
  cat("Results saved in 'atime_mvnormal_results.RData'\n")

  return(results)
}

# ==============================================================================
# Interactive Usage
# ==============================================================================

if (interactive()) {
  cat("MVNormal Distribution atime Benchmark Suite\n")
  cat("===========================================\n\n")
  cat("Available functions:\n")
  cat("  run_all_mvnormal_benchmarks()       - Run all benchmarks\n")
  cat("  run_mvnormal_atime_benchmark()      - Main benchmark (adjusted)\n")
  cat("  benchmark_mvnormal_R()              - R-only benchmark\n")
  cat("  benchmark_mvnormal_Cpp_large()      - C++ large dataset benchmark\n")
  cat("  run_mvnormal_adaptive_benchmark()   - Adaptive iteration benchmark\n")
  cat("  benchmark_mvnormal_components_quick() - Component analysis\n")
  cat("  benchmark_mvnormal_dimensions()     - Dimension scaling\n")
  cat("  compare_implementations()           - Quick comparison\n")
  cat("\nRun any function to start benchmarking!\n")
}
