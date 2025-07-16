#' Benchmark Integration Functions for atime Framework
#'
#' This file provides adapter functions for integrating Dirichlet Process
#' models with the atime benchmarking framework.
#'
#' @name benchmark_integration

#' Create standardized prior parameters for benchmarking
#'
#' @param dimensions Number of dimensions
#' @param model_name Covariance model name
#' @return List of prior parameters
#' @export
prepare_benchmark_parameters <- function(dimensions, model_name) {
  # Handle univariate models
  if (model_name %in% c("E", "V") && dimensions > 1) {
    stop("Models E and V are only for univariate data (dimensions = 1)")
  }
  
  # Handle multivariate-only models
  if (model_name %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI") && dimensions == 1) {
    stop(paste("Model", model_name, "is only for multivariate data (dimensions > 1)"))
  }
  
  # Base parameters
  prior_params <- list(
    mu0 = rep(0, dimensions),
    kappa0 = 1,
    nu = dimensions + 1,
    Lambda = diag(dimensions),
    covModel = model_name
  )
  
  return(prior_params)
}

#' Generate benchmark data for testing
#'
#' @param n Number of samples
#' @param d Number of dimensions
#' @param seed Random seed for reproducibility
#' @return Matrix of benchmark data
#' @export
generate_benchmark_data <- function(n, d, seed = 42) {
  set.seed(seed)
  
  if (d == 1) {
    # Univariate case - ensure it's a matrix
    data <- matrix(rnorm(n, mean = 0, sd = 1), ncol = 1)
  } else {
    # Multivariate case - create mixture of components
    k_clusters <- min(3, n %/% 5)  # Adaptive number of clusters
    cluster_sizes <- rep(n %/% k_clusters, k_clusters)
    cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)
    
    # Well-separated cluster means
    data_points <- c()
    for (i in 1:k_clusters) {
      mean_vec <- rep(0, d)
      mean_vec[1] <- (i - 2) * 2  # Separate along first dimension
      if (d > 1) mean_vec[2] <- (i - 2) * 1.5  # Separate along second dimension
      
      # Common covariance matrix
      sigma <- diag(d) * 0.5
      
      # Generate cluster data
      cluster_data <- mvtnorm::rmvnorm(cluster_sizes[i], mean_vec, sigma)
      data_points <- rbind(data_points, cluster_data)
    }
    
    data <- data_points[1:n, ]
  }
  
  return(data)
}

#' Collect standardized performance metrics from benchmark results
#'
#' @param result DP object result from benchmark
#' @param execution_time Execution time (optional)
#' @param memory_usage Memory usage (optional)
#' @return List of standardized metrics
#' @export
collect_benchmark_metrics <- function(result, execution_time = NULL, memory_usage = NULL) {
  # Extract standardized metrics from DP object
  metrics <- list(
    n_clusters = result$numberClusters,
    n_data_points = length(result$clusterLabels),
    dimensions = if (is.matrix(result$data)) ncol(result$data) else 1,
    model_type = if (!is.null(result$mixingDistribution$priorParameters$covModel)) {
      result$mixingDistribution$priorParameters$covModel
    } else "unknown",
    iterations = if (!is.null(result$iterations)) result$iterations else 0
  )
  
  # Add timing if available
  if (!is.null(execution_time)) {
    metrics$execution_time <- execution_time
  }
  
  # Add memory if available
  if (!is.null(memory_usage)) {
    metrics$memory_usage <- memory_usage
  }
  
  return(metrics)
}

#' Run optimized DP benchmark for atime integration
#'
#' @param data_matrix Input data matrix
#' @param model_name Covariance model name
#' @param mcmc_iterations Number of MCMC iterations (default: 10 for speed)
#' @param initial_clusters Number of initial clusters (default: 1)
#' @param return_format What to return ("metrics", "object", "clusters")
#' @return Benchmark result in specified format
#' @export
run_dp_benchmark <- function(data_matrix, model_name, mcmc_iterations = 10, 
                             initial_clusters = 1, return_format = "metrics") {
  
  # Validate inputs
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  
  dimensions <- ncol(data_matrix)
  
  # Create prior parameters
  prior_params <- prepare_benchmark_parameters(dimensions, model_name)
  
  # Create and run DP
  md <- MvnormalCreate(prior_params)
  dp <- DirichletProcessCreate(data_matrix, md)
  dp <- Initialise(dp, numInitialClusters = initial_clusters)
  
  # Run MCMC
  result <- Fit(dp, mcmc_iterations, progressBar = FALSE)
  
  # Return in requested format
  switch(return_format,
    "metrics" = {
      # Return simple metrics for atime compatibility
      c(
        n_clusters = result$numberClusters,
        n_points = length(result$clusterLabels),
        dimensions = dimensions
      )
    },
    "object" = {
      # Return full object
      result
    },
    "clusters" = {
      # Return just cluster count (simplest)
      result$numberClusters
    },
    {
      # Default: return cluster count
      result$numberClusters
    }
  )
}

#' Optimized atime benchmark for covariance models
#'
#' @param max_n Maximum sample size to test
#' @param dimensions Number of dimensions to use
#' @param models Vector of model names to benchmark
#' @param mcmc_iter Number of MCMC iterations (keep low for speed)
#' @param repetitions Number of benchmark repetitions
#' @return atime benchmark results
#' @export
run_optimized_atime_benchmark <- function(max_n = 100, dimensions = 2, 
                                          models = c("FULL", "EII"), 
                                          mcmc_iter = 5, repetitions = 3) {
  
  # Load atime
  if (!requireNamespace("atime", quietly = TRUE)) {
    stop("atime package is required for benchmarking")
  }
  
  # Filter models based on dimensions
  valid_models <- models
  if (dimensions == 1) {
    # Only E, V, and FULL work with univariate data
    valid_models <- intersect(models, c("E", "V", "FULL"))
  } else {
    # Remove E, V for multivariate data
    valid_models <- setdiff(models, c("E", "V"))
  }
  
  if (length(valid_models) == 0) {
    stop("No valid models for the specified dimensions")
  }
  
  # Create expressions dynamically
  expressions <- list()
  
  # Setup expression
  expressions$setup <- substitute({
    data_subset <- generate_benchmark_data(N, dimensions, seed = 42)
  }, list(dimensions = dimensions))
  
  # Create benchmark expressions for each model
  for (model in valid_models) {
    expressions[[model]] <- substitute({
      run_dp_benchmark(data_subset, model_name, mcmc_iterations = mcmc_iter, 
                       return_format = "clusters")
    }, list(model_name = model, mcmc_iter = mcmc_iter))
  }
  
  # Add times parameter
  expressions$times <- repetitions
  
  # Create N sequence
  n_sequence <- unique(c(10, 20, min(50, max_n), max_n))
  expressions$N <- n_sequence
  
  # Run benchmark
  result <- do.call(atime::atime, expressions)
  
  return(result)
}

#' Quick benchmark test for development
#'
#' @return atime benchmark results for testing
#' @export
quick_benchmark_test <- function() {
  run_optimized_atime_benchmark(
    max_n = 50, 
    dimensions = 2, 
    models = c("FULL", "EII"), 
    mcmc_iter = 3, 
    repetitions = 2
  )
}