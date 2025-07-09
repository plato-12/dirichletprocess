# debug_cluster_differences.R
# Script to debug cluster count differences between R and C++ implementations

library(dirichletprocess)
library(dplyr)

# Function to generate reproducible Gaussian mixture data
generate_gaussian_mixture <- function(n, k_true = 3, seed = 123) {
  set.seed(seed)

  # Generate data from k_true well-separated Gaussian clusters
  n_per_cluster <- n / k_true
  data <- numeric(n)

  # Create well-separated clusters
  means <- seq(-5, 5, length.out = k_true)
  sd <- 0.5

  for (i in 1:k_true) {
    start_idx <- floor((i-1) * n_per_cluster) + 1
    end_idx <- floor(i * n_per_cluster)
    data[start_idx:end_idx] <- rnorm(end_idx - start_idx + 1, means[i], sd)
  }

  # Shuffle to remove ordering
  data <- sample(data)

  return(list(data = data, true_k = k_true, true_means = means))
}

# Function to run single comparison
debug_single_comparison <- function(n_obs = 100, n_iter = 50, k_true = 3, seed = 42) {
  cat("\n=== Debugging Cluster Differences ===\n")
  cat(sprintf("Data: n=%d, true_k=%d\n", n_obs, k_true))
  cat(sprintf("MCMC: iterations=%d\n", n_iter))

  # Generate test data
  data_info <- generate_gaussian_mixture(n_obs, k_true, seed)
  data <- data_info$data

  # 1. Run R implementation
  cat("\n--- R Implementation ---\n")
  set.seed(100)  # Set seed for MCMC
  set_use_cpp(FALSE)

  dp_r <- DirichletProcessGaussian(data)
  cat("Initial clusters (R):", dp_r$numberClusters, "\n")
  cat("Initial alpha (R):", dp_r$alpha, "\n")

  # Run fit
  dp_r_fit <- Fit(dp_r, n_iter, progressBar = FALSE)

  cat("Final clusters (R):", dp_r_fit$numberClusters, "\n")
  cat("Final alpha (R):", dp_r_fit$alpha, "\n")
  cat("Cluster sizes (R):", sort(table(dp_r_fit$clusterLabels)), "\n")

  # 2. Run C++ implementation
  cat("\n--- C++ Implementation ---\n")
  set.seed(100)  # Same seed
  set_use_cpp(TRUE)

  dp_cpp <- DirichletProcessGaussian(data)
  cat("Initial clusters (C++):", dp_cpp$numberClusters, "\n")
  cat("Initial alpha (C++):", dp_cpp$alpha, "\n")

  # Run fit
  dp_cpp_fit <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  cat("Final clusters (C++):", dp_cpp_fit$numberClusters, "\n")
  cat("Final alpha (C++):", dp_cpp_fit$alpha, "\n")
  cat("Cluster sizes (C++):", sort(table(dp_cpp_fit$clusterLabels)), "\n")

  # 3. Compare cluster evolution
  cat("\n--- Cluster Evolution ---\n")

  # Extract cluster counts over iterations
  r_cluster_evolution <- sapply(dp_r_fit$labelsChain, function(x) length(unique(x)))
  cpp_cluster_evolution <- sapply(dp_cpp_fit$labelsChain, function(x) length(unique(x)))

  cat("R cluster evolution (first 10):", head(r_cluster_evolution, 10), "\n")
  cat("C++ cluster evolution (first 10):", head(cpp_cluster_evolution, 10), "\n")

  # 4. Compare likelihood evolution
  cat("\n--- Likelihood Evolution ---\n")
  cat("R likelihood (last 5):", tail(dp_r_fit$likelihoodChain, 5), "\n")
  cat("C++ likelihood (last 5):", tail(dp_cpp_fit$likelihoodChain, 5), "\n")

  # 5. Check hyperparameters
  cat("\n--- Hyperparameters ---\n")
  cat("R mixing distribution class:", class(dp_r_fit$mixingDistribution), "\n")
  cat("C++ mixing distribution class:", class(dp_cpp_fit$mixingDistribution), "\n")

  if (!is.null(dp_r_fit$mixingDistribution$priors)) {
    cat("R priors:", paste(names(dp_r_fit$mixingDistribution$priors), "=",
                           unlist(dp_r_fit$mixingDistribution$priors), collapse=", "), "\n")
  }
  if (!is.null(dp_cpp_fit$mixingDistribution$priors)) {
    cat("C++ priors:", paste(names(dp_cpp_fit$mixingDistribution$priors), "=",
                             unlist(dp_cpp_fit$mixingDistribution$priors), collapse=", "), "\n")
  }

  # Return results for further analysis
  return(list(
    r = dp_r_fit,
    cpp = dp_cpp_fit,
    data = data,
    true_k = k_true
  ))
}

# Function to check algorithm parameters
check_algorithm_parameters <- function() {
  cat("\n=== Checking Algorithm Parameters ===\n")

  # Create a simple dataset
  data <- rnorm(50)

  # Check R implementation settings
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(data)

  cat("R Implementation:\n")
  cat("  m_auxiliary (if exists):",
      ifelse(exists("m_auxiliary", dp_r), dp_r$m_auxiliary, "Not found"), "\n")
  cat("  mhDraws:", dp_r$mhDraws, "\n")
  cat("  alpha:", dp_r$alpha, "\n")

  # Check C++ implementation settings
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(data)

  cat("\nC++ Implementation:\n")
  cat("  m_auxiliary (if exists):",
      ifelse(exists("m_auxiliary", dp_cpp), dp_cpp$m_auxiliary, "Not found"), "\n")
  cat("  mhDraws:", dp_cpp$mhDraws, "\n")
  cat("  alpha:", dp_cpp$alpha, "\n")

  # Check C++ status
  cat("\nC++ Status:\n")
  cpp_status <- get_cpp_status()
  print(cpp_status)
}

# Function to test with controlled initialization
test_controlled_initialization <- function(n_obs = 100, n_iter = 50) {
  cat("\n=== Testing Controlled Initialization ===\n")

  # Generate simple data
  data <- c(rnorm(50, -3, 0.3), rnorm(50, 3, 0.3))

  # Test with fixed initial number of clusters
  for (init_clusters in c(1, 3, 5)) {
    cat(sprintf("\n--- Initial clusters: %d ---\n", init_clusters))

    # R implementation
    set.seed(200)
    set_use_cpp(FALSE)
    dp_r <- DirichletProcessGaussian(data)
    dp_r <- Initialise(dp_r, numInitialClusters = init_clusters)
    dp_r_fit <- Fit(dp_r, n_iter, progressBar = FALSE)

    # C++ implementation
    set.seed(200)
    set_use_cpp(TRUE)
    dp_cpp <- DirichletProcessGaussian(data)
    dp_cpp <- Initialise(dp_cpp, numInitialClusters = init_clusters)
    dp_cpp_fit <- Fit(dp_cpp, n_iter, progressBar = FALSE)

    cat(sprintf("R: %d -> %d clusters\n", init_clusters, dp_r_fit$numberClusters))
    cat(sprintf("C++: %d -> %d clusters\n", init_clusters, dp_cpp_fit$numberClusters))
  }
}

# Main debugging execution
if (interactive()) {
  cat("Starting cluster difference debugging...\n")

  # 1. Check algorithm parameters
  check_algorithm_parameters()

  # 2. Run detailed comparison
  results <- debug_single_comparison(n_obs = 200, n_iter = 100, k_true = 3)

  # 3. Test with controlled initialization
  test_controlled_initialization()

  # 4. Test multiple scenarios
  cat("\n=== Testing Multiple Scenarios ===\n")
  scenarios <- expand.grid(
    n_obs = c(100, 200),
    k_true = c(2, 3),
    n_iter = c(50, 100)
  )

  scenario_results <- list()
  for (i in 1:nrow(scenarios)) {
    cat(sprintf("\nScenario %d/%d\n", i, nrow(scenarios)))
    s <- scenarios[i,]

    res <- debug_single_comparison(s$n_obs, s$n_iter, s$k_true, seed = i * 100)

    scenario_results[[i]] <- data.frame(
      scenario = i,
      n_obs = s$n_obs,
      k_true = s$k_true,
      n_iter = s$n_iter,
      r_clusters = res$r$numberClusters,
      cpp_clusters = res$cpp$numberClusters,
      cluster_diff = abs(res$r$numberClusters - res$cpp$numberClusters)
    )
  }

  # Summary
  summary_df <- bind_rows(scenario_results)
  cat("\n=== Summary of All Scenarios ===\n")
  print(summary_df)

  cat("\nAverage cluster difference:", mean(summary_df$cluster_diff), "\n")
  cat("Max cluster difference:", max(summary_df$cluster_diff), "\n")
}
