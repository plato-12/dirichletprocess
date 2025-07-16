# Optimized Covariance Models Benchmark
# =====================================
# 
# This script provides an optimized version of the comprehensive covariance models
# benchmark that works efficiently with the atime framework.
#
# Key improvements:
# - Reduced MCMC iterations for faster execution
# - Standardized return formats for atime compatibility
# - Proper error handling and validation
# - Support for both univariate (E, V) and multivariate models
#
# Usage:
#   source("benchmark/atime/benchmark-covariance-models-optimized.R")
#   result <- run_fast_covariance_benchmark()

# Required libraries
library(dirichletprocess)
if (!require(atime, quietly = TRUE)) {
  stop("atime package is required. Install with: install.packages('atime')")
}
if (!require(mvtnorm, quietly = TRUE)) {
  stop("mvtnorm package is required. Install with: install.packages('mvtnorm')")
}

# Source benchmark integration functions
if (!exists("prepare_benchmark_parameters")) {
  source("R/benchmark_integration.R")
}

# ==========================================
# OPTIMIZED BENCHMARK CONFIGURATIONS
# ==========================================

# Fast benchmark configuration
FAST_CONFIG <- list(
  mcmc_iterations = 5,     # Very fast for development/testing
  repetitions = 3,
  max_samples = 50
)

# Standard benchmark configuration  
STANDARD_CONFIG <- list(
  mcmc_iterations = 20,    # Reasonable for CI/regular testing
  repetitions = 5,
  max_samples = 100
)

# Comprehensive benchmark configuration
COMPREHENSIVE_CONFIG <- list(
  mcmc_iterations = 50,    # More thorough for research/publication
  repetitions = 10,
  max_samples = 200
)

# ==========================================
# MAIN BENCHMARK FUNCTIONS
# ==========================================

#' Run fast covariance model benchmark
#'
#' @param config Configuration list (FAST_CONFIG, STANDARD_CONFIG, COMPREHENSIVE_CONFIG)
#' @param dimensions Number of dimensions to test
#' @param models Vector of models to test (NULL for auto-selection)
#' @return atime benchmark results
#' @export
run_fast_covariance_benchmark <- function(config = FAST_CONFIG, dimensions = 2, models = NULL) {
  
  cat("=== Optimized Covariance Models Benchmark ===\n")
  cat("Configuration:", deparse(substitute(config)), "\n")
  cat("Dimensions:", dimensions, "\n")
  cat("MCMC iterations:", config$mcmc_iterations, "\n")
  cat("Repetitions:", config$repetitions, "\n")
  cat("Max samples:", config$max_samples, "\n\n")
  
  # Auto-select models based on dimensions
  if (is.null(models)) {
    if (dimensions == 1) {
      models <- c("E", "V", "FULL")
    } else {
      models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
    }
  }
  
  cat("Models to test:", paste(models, collapse = ", "), "\n\n")
  
  # Run optimized benchmark
  result <- run_optimized_atime_benchmark(
    max_n = config$max_samples,
    dimensions = dimensions,
    models = models,
    mcmc_iter = config$mcmc_iterations,
    repetitions = config$repetitions
  )
  
  cat("=== Benchmark Complete ===\n")
  
  return(result)
}

#' Run univariate models benchmark
#'
#' @param config Configuration list
#' @return atime benchmark results for E, V, and FULL models
#' @export  
run_univariate_benchmark <- function(config = FAST_CONFIG) {
  cat("=== Univariate Models Benchmark (E, V, FULL) ===\n")
  
  result <- run_fast_covariance_benchmark(
    config = config,
    dimensions = 1,
    models = c("E", "V", "FULL")
  )
  
  return(result)
}

#' Run multivariate models benchmark
#'
#' @param config Configuration list
#' @param dimensions Number of dimensions (must be > 1)
#' @return atime benchmark results for multivariate models
#' @export
run_multivariate_benchmark <- function(config = STANDARD_CONFIG, dimensions = 2) {
  if (dimensions <= 1) {
    stop("Multivariate benchmark requires dimensions > 1")
  }
  
  cat("=== Multivariate Models Benchmark ===\n")
  
  result <- run_fast_covariance_benchmark(
    config = config,
    dimensions = dimensions,
    models = c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
  )
  
  return(result)
}

#' Run comprehensive benchmark across multiple dimensions
#'
#' @param config Configuration list
#' @param dimensions_list Vector of dimensions to test
#' @return List of benchmark results
#' @export
run_comprehensive_benchmark <- function(config = STANDARD_CONFIG, 
                                       dimensions_list = c(1, 2, 5)) {
  
  cat("=== Comprehensive Multi-Dimensional Benchmark ===\n")
  
  results <- list()
  
  for (d in dimensions_list) {
    cat(sprintf("\n--- Testing %dD models ---\n", d))
    
    if (d == 1) {
      result <- run_univariate_benchmark(config)
      results[[paste0("d", d)]] <- result
    } else {
      result <- run_multivariate_benchmark(config, d)
      results[[paste0("d", d)]] <- result
    }
  }
  
  cat("\n=== Comprehensive Benchmark Complete ===\n")
  
  return(results)
}

# ==========================================
# ANALYSIS AND REPORTING FUNCTIONS
# ==========================================

#' Analyze benchmark results
#'
#' @param benchmark_result atime benchmark result object
#' @return Summary analysis
#' @export
analyze_benchmark_results <- function(benchmark_result) {
  measurements <- benchmark_result$measurements
  
  # Calculate summary statistics
  summary_stats <- measurements[, .(
    mean_time = mean(median),
    min_time = min(median),
    max_time = max(median),
    std_time = sd(median),
    mean_memory = mean(kilobytes, na.rm = TRUE)
  ), by = expr.name]
  
  # Find fastest and slowest models
  fastest_model <- summary_stats[which.min(mean_time)]$expr.name
  slowest_model <- summary_stats[which.max(mean_time)]$expr.name
  
  # Calculate performance ratios
  baseline_time <- summary_stats[expr.name == "FULL"]$mean_time
  if (length(baseline_time) > 0) {
    summary_stats[, performance_ratio := mean_time / baseline_time]
  }
  
  analysis <- list(
    summary_stats = summary_stats,
    fastest_model = fastest_model,
    slowest_model = slowest_model,
    total_models_tested = nrow(summary_stats),
    sample_sizes_tested = unique(measurements$N)
  )
  
  return(analysis)
}

#' Print benchmark summary
#'
#' @param benchmark_result atime benchmark result object
#' @export
print_benchmark_summary <- function(benchmark_result) {
  analysis <- analyze_benchmark_results(benchmark_result)
  
  cat("=== Benchmark Summary ===\n")
  cat("Models tested:", analysis$total_models_tested, "\n")
  cat("Sample sizes:", paste(analysis$sample_sizes_tested, collapse = ", "), "\n")
  cat("Fastest model:", analysis$fastest_model, "\n")
  cat("Slowest model:", analysis$slowest_model, "\n\n")
  
  cat("Performance Summary:\n")
  print(analysis$summary_stats[order(mean_time)])
  
  return(analysis)
}

# ==========================================
# CONVENIENCE FUNCTIONS
# ==========================================

#' Quick test of all working models
#'
#' @export
quick_test_all_models <- function() {
  cat("=== Quick Test of All Models ===\n")
  
  # Test univariate
  cat("Testing univariate models...\n")
  univariate_result <- run_univariate_benchmark(FAST_CONFIG)
  
  # Test multivariate
  cat("Testing multivariate models...\n") 
  multivariate_result <- run_multivariate_benchmark(FAST_CONFIG, dimensions = 2)
  
  # Print summaries
  cat("\nUnivariate Results:\n")
  print_benchmark_summary(univariate_result)
  
  cat("\nMultivariate Results:\n")
  print_benchmark_summary(multivariate_result)
  
  return(list(
    univariate = univariate_result,
    multivariate = multivariate_result
  ))
}

# ==========================================
# EXAMPLE USAGE
# ==========================================

if (FALSE) {
  # Example usage (set to TRUE to run)
  
  # Quick test
  quick_results <- quick_test_all_models()
  
  # Standard benchmark
  standard_result <- run_fast_covariance_benchmark(STANDARD_CONFIG, dimensions = 2)
  print_benchmark_summary(standard_result)
  
  # Comprehensive multi-dimensional
  comprehensive_results <- run_comprehensive_benchmark(STANDARD_CONFIG, c(1, 2, 5))
  
  # Individual model types
  univariate_only <- run_univariate_benchmark(COMPREHENSIVE_CONFIG)
  multivariate_only <- run_multivariate_benchmark(COMPREHENSIVE_CONFIG, dimensions = 5)
}

cat("Optimized benchmark functions loaded successfully.\n")
cat("Try: quick_test_all_models()\n")