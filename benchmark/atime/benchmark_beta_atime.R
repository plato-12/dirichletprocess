# benchmark_beta_atime.R - FIXED VERSION
# Benchmark Beta Distribution using atime package
# Based on Neal (2000) and Escobar & West (1995) algorithms

library(dirichletprocess)
library(atime)
library(ggplot2)

# Fix 1: Update the Likelihood.beta function to handle both old and new theta formats
fix_beta_likelihood <- function() {
  assignInNamespace("Likelihood.beta", function(mdObj, x, theta) {
    maxT <- mdObj$maxT
    x <- as.vector(x, "numeric")

    # Handle both indexed and named theta formats
    if (is.list(theta) && !is.null(names(theta)) && "mu" %in% names(theta)) {
      # New format with named components
      mu <- theta$mu
      nu <- theta$nu
    } else {
      # Old format with indexed components
      mu <- as.numeric(theta[[1]][, , , drop = TRUE])
      tau <- as.numeric(theta[[2]][, , , drop = TRUE])
      # Convert tau to nu for consistency
      nu <- tau
    }

    # Ensure we have values
    if (length(mu) == 0 || length(nu) == 0) {
      return(numeric(length(x)))
    }

    # Ensure mu and nu are numeric vectors
    mu <- as.numeric(mu)
    nu <- as.numeric(nu)

    # Recycle parameters if needed
    mu <- rep_len(mu, length(x))
    nu <- rep_len(nu, length(x))

    # Calculate likelihood
    y <- numeric(length(x))
    for (i in seq_along(x)) {
      # Validate parameters
      if (is.na(mu[i]) || is.na(nu[i]) || mu[i] <= 0 || mu[i] >= maxT || nu[i] <= 0) {
        y[i] <- 1e-300
        next
      }

      a <- (mu[i] * nu[i]) / maxT
      b <- (1 - mu[i]/maxT) * nu[i]

      # Ensure valid beta parameters
      if (a <= 0 || b <= 0 || !is.finite(a) || !is.finite(b)) {
        y[i] <- 1e-300
        next
      }

      # Calculate likelihood
      if (x[i] >= 0 && x[i] <= maxT) {
        y[i] <- (1/maxT) * dbeta(x[i]/maxT, a, b)
      } else {
        y[i] <- 1e-300
      }
    }

    return(as.numeric(y))
  }, ns = "dirichletprocess")
}

# Fix 2: Update Initialise.beta to use correct theta structure
fix_beta_initialization <- function() {
  assignInNamespace("Initialise.beta", function(dpObj, posterior = TRUE, verbose = TRUE, ...) {

    # Ensure all points start in cluster 1
    dpObj$clusterLabels <- rep(1, dpObj$n)
    dpObj$numberClusters <- 1
    dpObj$pointsPerCluster <- numeric(dpObj$n)
    dpObj$pointsPerCluster[1] <- dpObj$n

    # Initialize parameters with correct array structure
    if (posterior) {
      cluster_data <- matrix(dpObj$data, ncol = 1)
      post_draws <- PosteriorDraw(dpObj$mixingDistribution, cluster_data, n = 1)

      # Ensure parameters are in array format
      dpObj$clusterParameters <- list(
        array(as.numeric(post_draws$mu), dim = c(1, 1, 1)),
        array(as.numeric(post_draws$nu), dim = c(1, 1, 1))
      )
    } else {
      prior_draws <- PriorDraw(dpObj$mixingDistribution, 1)

      # Ensure parameters are in array format
      dpObj$clusterParameters <- list(
        array(as.numeric(prior_draws$mu), dim = c(1, 1, 1)),
        array(as.numeric(prior_draws$nu), dim = c(1, 1, 1))
      )
    }

    # Initialize auxiliary parameters for non-conjugate
    dpObj$m <- 3
    dpObj$aux <- vector("list", dpObj$m)
    for (j in seq_len(dpObj$m)) {
      aux_params <- PriorDraw(dpObj$mixingDistribution, 1)
      dpObj$aux[[j]] <- list(
        mu = as.numeric(aux_params$mu),
        nu = as.numeric(aux_params$nu)
      )
    }

    if (verbose) {
      cat("Initialised Dirichlet process with 1 cluster\n")
    }

    return(dpObj)
  }, ns = "dirichletprocess")
}

# Apply both fixes
fix_beta_likelihood()
fix_beta_initialization()

# ==============================================================================
# Setup Functions
# ==============================================================================

#' Generate synthetic Beta mixture data for benchmarking
#' Creates mixture of Beta distributions with known clusters
generate_beta_test_data <- function(n, k_clusters = 2, seed = 42) {
  set.seed(seed)

  # Equal sized clusters
  cluster_sizes <- rep(n %/% k_clusters, k_clusters)
  cluster_sizes[k_clusters] <- cluster_sizes[k_clusters] + (n %% k_clusters)

  # Well-separated Beta parameters (alpha, beta)
  # These create distinct shapes: left-skewed, right-skewed, symmetric
  shapes <- list(
    c(2, 8),   # Left-skewed
    c(8, 2),   # Right-skewed
    c(5, 5)    # Symmetric
  )[1:k_clusters]

  data <- numeric(n)
  idx <- 1

  for (i in 1:k_clusters) {
    cluster_data <- rbeta(cluster_sizes[i],
                          shape1 = shapes[[i]][1],
                          shape2 = shapes[[i]][2])
    data[idx:(idx + cluster_sizes[i] - 1)] <- cluster_data
    idx <- idx + cluster_sizes[i]
  }

  return(data)
}

# ==============================================================================
# Main atime Benchmark
# ==============================================================================

#' Run atime benchmark comparing R and C++ implementations for Beta
run_beta_atime_benchmark <- function() {

  cat("========================================\n")
  cat("Beta Distribution DP Benchmark (atime)\n")
  cat("========================================\n\n")

  # Check C++ availability for Beta
  cpp_available <- exists("_dirichletprocess_run_mcmc_cpp") &&
    exists("can_use_cpp")

  if (!cpp_available) {
    cat("⚠️  C++ implementation not available for Beta\n")
    cat("Only R implementation will be benchmarked\n\n")
  } else {
    # Test if C++ supports Beta
    test_data <- rbeta(10, 2, 2)
    test_dp <- DirichletProcessBeta(test_data, verbose = FALSE)
    cpp_supports_beta <- tryCatch({
      can_use_cpp(test_dp)
    }, error = function(e) FALSE)

    if (!cpp_supports_beta) {
      cat("⚠️  C++ implementation does not support Beta distribution yet\n")
      cat("Only R implementation will be benchmarked\n\n")
      cpp_available <- FALSE
    } else {
      cat("✓ C++ implementation available for Beta\n\n")
    }
  }

  # Define expression list for atime
  expr_list <- list()

  # R implementation expression
  expr_list$R_implementation <- quote({
    set_use_cpp(FALSE)
    set.seed(123)
    dp <- DirichletProcessBeta(data, verbose = FALSE)
    dp <- Fit(dp, iterations, progressBar = FALSE)
  })

  # C++ implementation expression (if available)
  if (cpp_available) {
    expr_list$Cpp_implementation <- quote({
      set_use_cpp(TRUE)
      set.seed(123)
      dp <- DirichletProcessBeta(data, verbose = FALSE)
      dp <- Fit(dp, iterations, progressBar = FALSE)
    })
  }

  # Run atime benchmark
  atime_result <- atime::atime(
    N = as.integer(10^seq(1.5, 3.5, by = 0.25)),  # 30 to 3000+ observations
    setup = {
      data <- generate_beta_test_data(N, k_clusters = 2)
      iterations <- 100  # Fixed number of iterations
    },
    expr.list = expr_list,
    seconds.limit = 10,  # Stop if any expression takes more than 10 seconds
    verbose = TRUE
  )

  # Print summary
  print(atime_result)

  # Create and save plot
  p <- plot(atime_result)
  print(p)

  # Save results
  cat("\nSaving results...\n")
  save(atime_result, file = "atime_beta_results.RData")
  ggsave("atime_beta_benchmark.png", plot = p, width = 10, height = 8, dpi = 300)

  return(atime_result)
}

# ==============================================================================
# Component-level Benchmarking with atime
# ==============================================================================

#' Benchmark individual components for Beta using atime
benchmark_beta_components <- function() {

  cat("\n\n==========================================\n")
  cat("Beta Component-level Benchmarking (atime)\n")
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
        mdObj <- BetaMixtureCreate(c(2, 8), mhStepSize = c(0.1, 0.1))
        # Create theta with proper array structure
        theta <- list(
          array(0.5, dim = c(1, 1, 1)),
          array(10, dim = c(1, 1, 1))
        )
        # Calculate likelihood for all data points
        for (i in 1:10) {
          lik <- Likelihood(mdObj, data, theta)
        }
      })
    } else if (comp == "PriorDraw") {
      expr_list_comp[["Prior_sampling"]] <- quote({
        mdObj <- BetaMixtureCreate(c(2, 8), mhStepSize = c(0.1, 0.1))
        # Draw multiple prior samples
        for (i in 1:10) {
          prior_samples <- PriorDraw(mdObj, 10)
        }
      })
    } else if (comp == "PosteriorDraw") {
      expr_list_comp[["Posterior_sampling"]] <- quote({
        mdObj <- BetaMixtureCreate(c(2, 8), mhStepSize = c(0.1, 0.1))
        # Draw posterior samples using MH
        post_samples <- PosteriorDraw(mdObj, matrix(data, ncol = 1), n = 10)
      })
    }

    # Run component benchmark
    comp_result <- atime::atime(
      N = as.integer(10^seq(2, 3.5, by = 0.5)),
      setup = {
        data <- generate_beta_test_data(N)
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

#' Analyze memory scaling for Beta distribution
analyze_beta_memory_scaling <- function() {

  cat("\n\n========================================\n")
  cat("Beta Memory Scaling Analysis (atime)\n")
  cat("========================================\n\n")

  # Focus on memory measurement
  expr_list_mem <- list()

  expr_list_mem$R_memory <- quote({
    set_use_cpp(FALSE)
    dp <- DirichletProcessBeta(data, verbose = FALSE)
    dp <- Fit(dp, 50, progressBar = FALSE)
    # Force garbage collection to get accurate memory usage
    gc()
  })

  # Run memory-focused benchmark
  memory_result <- atime::atime(
    N = as.integer(10^seq(2, 4, by = 0.5)),  # Up to 10,000 observations
    setup = {
      data <- generate_beta_test_data(N)
    },
    expr.list = expr_list_mem,
    seconds.limit = 20
  )

  print(memory_result)
  plot(memory_result)

  return(memory_result)
}

# ==============================================================================
# Main Execution with Fixed Initialization
# ==============================================================================

if (interactive()) {
  cat("Beta Distribution atime Benchmark Suite\n")
  cat("======================================\n\n")
  cat("Available benchmarks:\n")
  cat("  run_beta_atime_benchmark()      - Main scaling comparison\n")
  cat("  benchmark_beta_components()     - Component-level analysis\n")
  cat("  analyze_beta_memory_scaling()   - Memory usage scaling\n")
  cat("\nRun any function to start benchmarking!\n")

  # Quick test to verify setup
  cat("\nRunning quick verification test...\n")
  test_data <- generate_beta_test_data(100)

  set_use_cpp(FALSE)
  dp_test <- DirichletProcessBeta(test_data, verbose = FALSE)
  dp_test <- Fit(dp_test, 10, progressBar = FALSE)
  cat(sprintf("✓ R implementation works (found %d clusters)\n", dp_test$numberClusters))

  cat("\nReady for benchmarking!\n")
}
