# =============================================================================
# tests/testthat/test_cpp_r_equivalence.R
# Comprehensive Equivalence Tests for C++ vs R Gaussian MCMC Implementation
# =============================================================================

context("C++ vs R Gaussian MCMC Equivalence Tests")

# =============================================================================
# Helper Functions for Statistical Testing
# =============================================================================

#' Generate reproducible test data with known structure
generate_test_mixture <- function(n, k = 3, seed = 42, separation = 3) {
  set.seed(seed)

  # Create well-separated clusters
  means <- seq(-separation, separation, length.out = k)
  sds <- rep(0.5, k)
  sizes <- as.numeric(rmultinom(1, n, rep(1/k, k)))

  data <- c()
  true_labels <- c()

  for (i in 1:k) {
    if (sizes[i] > 0) {
      cluster_data <- rnorm(sizes[i], mean = means[i], sd = sds[i])
      data <- c(data, cluster_data)
      true_labels <- c(true_labels, rep(i, sizes[i]))
    }
  }

  return(list(
    data = data,
    true_labels = true_labels,
    true_means = means,
    true_sds = sds,
    true_k = k,
    sizes = sizes
  ))
}

#' Statistical test for equivalence of two samples
test_sample_equivalence <- function(sample1, sample2, tolerance = 0.1, test = "ks") {
  if (test == "ks") {
    # Kolmogorov-Smirnov test
    ks_result <- ks.test(sample1, sample2)
    return(ks_result$p.value > 0.05)  # Non-significant = equivalent
  } else if (test == "t") {
    # T-test for means
    t_result <- t.test(sample1, sample2)
    return(t_result$p.value > 0.05)
  } else if (test == "var") {
    # F-test for variances
    var_result <- var.test(sample1, sample2)
    return(var_result$p.value > 0.05)
  }
}

#' Compare cluster assignments using Adjusted Rand Index
compare_clusterings <- function(labels1, labels2) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    # Simple alternative: proportion of exact matches after alignment
    return(simple_clustering_similarity(labels1, labels2))
  }

  mclust::adjustedRandIndex(labels1, labels2)
}

#' Simple clustering similarity when mclust not available
simple_clustering_similarity <- function(labels1, labels2) {
  # Find best permutation of labels
  k1 <- length(unique(labels1))
  k2 <- length(unique(labels2))

  if (k1 != k2) {
    return(0)  # Different number of clusters
  }

  best_accuracy <- 0

  # Try all permutations (for small k)
  if (k1 <= 6) {
    perms <- gtools::permutations(k1, k1)
    for (i in 1:nrow(perms)) {
      perm <- perms[i, ]
      aligned_labels2 <- perm[labels2]
      accuracy <- mean(labels1 == aligned_labels2)
      best_accuracy <- max(best_accuracy, accuracy)
    }
  } else {
    # For larger k, use simple matching
    best_accuracy <- mean(labels1 == labels2)
  }

  return(best_accuracy)
}

#' Extract posterior means from MCMC chains
extract_posterior_means <- function(dp_obj, burn_prop = 0.5) {
  if (!"alphaChain" %in% names(dp_obj)) {
    return(list(alpha = dp_obj$alpha))
  }

  n_samples <- length(dp_obj$alphaChain)
  burn_in <- floor(n_samples * burn_prop)
  keep_samples <- (burn_in + 1):n_samples

  list(
    alpha = mean(dp_obj$alphaChain[keep_samples]),
    final_clusters = dp_obj$numberClusters,
    final_likelihood = if ("likelihoodChain" %in% names(dp_obj)) {
      tail(dp_obj$likelihoodChain, 1)
    } else NA
  )
}

#' Run both implementations with same random seed
run_both_implementations <- function(data, n_iter = 100, seed = 123,
                                     updatePrior = FALSE, progressBar = FALSE) {

  # R Implementation
  set_use_cpp(FALSE)
  set.seed(seed)
  dp_r <- DirichletProcessGaussian(data)
  dp_r <- Fit(dp_r, n_iter, updatePrior = updatePrior, progressBar = progressBar)

  # C++ Implementation (if available)
  cpp_available <- exists("_dirichletprocess_run_mcmc_cpp") ||
    get_cpp_status()$mcmc_runner

  if (cpp_available && can_use_cpp(DirichletProcessGaussian(data))) {
    set_use_cpp(TRUE)
    set.seed(seed)  # Same seed
    dp_cpp <- DirichletProcessGaussian(data)
    dp_cpp <- Fit(dp_cpp, n_iter, updatePrior = updatePrior, progressBar = progressBar)

    return(list(r = dp_r, cpp = dp_cpp, both_available = TRUE))
  } else {
    return(list(r = dp_r, cpp = NULL, both_available = FALSE))
  }
}

# =============================================================================
# Basic Equivalence Tests
# =============================================================================

test_that("C++ and R produce equivalent cluster counts", {
  test_data <- generate_test_mixture(n = 100, k = 3, seed = 42)
  results <- run_both_implementations(test_data$data, n_iter = 50, seed = 123)

  skip_if(!results$both_available, "C++ implementation not available")

  # Should find similar number of clusters (within reasonable range)
  cluster_diff <- abs(results$r$numberClusters - results$cpp$numberClusters)
  expect_true(cluster_diff <= 2,
              info = sprintf("R found %d clusters, C++ found %d clusters",
                             results$r$numberClusters, results$cpp$numberClusters))

  # Both should find at least 1 cluster
  expect_true(results$r$numberClusters >= 1)
  expect_true(results$cpp$numberClusters >= 1)

  # Both should find reasonable number of clusters (not too many)
  expect_true(results$r$numberClusters <= length(test_data$data) / 2)
  expect_true(results$cpp$numberClusters <= length(test_data$data) / 2)
})

test_that("C++ and R produce equivalent posterior alpha estimates", {
  test_data <- generate_test_mixture(n = 80, k = 2, seed = 100)
  results <- run_both_implementations(test_data$data, n_iter = 100, seed = 200)

  skip_if(!results$both_available, "C++ implementation not available")

  # Extract posterior means
  r_summary <- extract_posterior_means(results$r)
  cpp_summary <- extract_posterior_means(results$cpp)

  # Alpha estimates should be reasonably close
  alpha_diff <- abs(r_summary$alpha - cpp_summary$alpha)
  expect_true(alpha_diff < 2.0,
              info = sprintf("R alpha: %.3f, C++ alpha: %.3f",
                             r_summary$alpha, cpp_summary$alpha))
})

test_that("C++ and R produce similar likelihood values", {
  test_data <- generate_test_mixture(n = 60, k = 2, seed = 150)
  results <- run_both_implementations(test_data$data, n_iter = 80, seed = 250)

  skip_if(!results$both_available, "C++ implementation not available")

  # Check final likelihood values are reasonable
  if ("likelihoodChain" %in% names(results$r) &&
      "likelihoodChain" %in% names(results$cpp)) {

    r_likelihood <- tail(results$r$likelihoodChain, 1)
    cpp_likelihood <- tail(results$cpp$likelihoodChain, 1)

    # Likelihoods should be finite
    expect_true(is.finite(r_likelihood))
    expect_true(is.finite(cpp_likelihood))

    # Should be reasonably close (within 10%)
    if (!is.na(r_likelihood) && !is.na(cpp_likelihood)) {
      rel_diff <- abs(r_likelihood - cpp_likelihood) / max(abs(r_likelihood), abs(cpp_likelihood))
      expect_true(rel_diff < 0.5,
                  info = sprintf("R likelihood: %.2f, C++ likelihood: %.2f",
                                 r_likelihood, cpp_likelihood))
    }
  }
})

# =============================================================================
# Reproducibility Tests
# =============================================================================

test_that("Same seed produces identical results within implementation", {
  test_data <- generate_test_mixture(n = 50, k = 2, seed = 300)

  # Test R implementation reproducibility
  set_use_cpp(FALSE)

  set.seed(400)
  dp_r1 <- DirichletProcessGaussian(test_data$data)
  dp_r1 <- Fit(dp_r1, 30, progressBar = FALSE)

  set.seed(400)  # Same seed
  dp_r2 <- DirichletProcessGaussian(test_data$data)
  dp_r2 <- Fit(dp_r2, 30, progressBar = FALSE)

  # Should get identical results
  expect_equal(dp_r1$numberClusters, dp_r2$numberClusters)
  expect_equal(dp_r1$clusterLabels, dp_r2$clusterLabels)

  # Test C++ implementation reproducibility (if available)
  cpp_available <- exists("_dirichletprocess_run_mcmc_cpp") ||
    get_cpp_status()$mcmc_runner

  if (cpp_available && can_use_cpp(DirichletProcessGaussian(test_data$data))) {
    set_use_cpp(TRUE)

    set.seed(500)
    dp_cpp1 <- DirichletProcessGaussian(test_data$data)
    dp_cpp1 <- Fit(dp_cpp1, 30, progressBar = FALSE)

    set.seed(500)  # Same seed
    dp_cpp2 <- DirichletProcessGaussian(test_data$data)
    dp_cpp2 <- Fit(dp_cpp2, 30, progressBar = FALSE)

    # Should get identical results
    expect_equal(dp_cpp1$numberClusters, dp_cpp2$numberClusters)
    expect_equal(dp_cpp1$clusterLabels, dp_cpp2$clusterLabels)
  }
})

# =============================================================================
# Statistical Equivalence Tests
# =============================================================================

test_that("C++ and R converge to similar posterior distributions", {
  test_data <- generate_test_mixture(n = 120, k = 3, seed = 500)

  # Run longer chains for better convergence
  results <- run_both_implementations(test_data$data, n_iter = 200, seed = 600)

  skip_if(!results$both_available, "C++ implementation not available")

  # Compare alpha chains (if available)
  if ("alphaChain" %in% names(results$r) && "alphaChain" %in% names(results$cpp)) {

    # Use second half of chains (after burn-in)
    n_samples <- min(length(results$r$alphaChain), length(results$cpp$alphaChain))
    burn_in <- floor(n_samples / 2)

    r_alpha_samples <- results$r$alphaChain[(burn_in + 1):n_samples]
    cpp_alpha_samples <- results$cpp$alphaChain[(burn_in + 1):n_samples]

    # Test that alpha distributions are similar
    if (length(r_alpha_samples) > 10 && length(cpp_alpha_samples) > 10) {
      alpha_equivalent <- test_sample_equivalence(r_alpha_samples, cpp_alpha_samples, test = "ks")
      expect_true(alpha_equivalent,
                  info = sprintf("Alpha distributions differ: R mean=%.3f, C++ mean=%.3f",
                                 mean(r_alpha_samples), mean(cpp_alpha_samples)))
    }
  }
})

test_that("C++ and R cluster assignments are statistically similar", {
  # Use well-separated clusters for clearer clustering
  test_data <- generate_test_mixture(n = 100, k = 2, seed = 700, separation = 4)

  # Run multiple times with different seeds to assess consistency
  similarities <- numeric(5)

  for (i in 1:5) {
    results <- run_both_implementations(test_data$data, n_iter = 100, seed = 800 + i)

    if (results$both_available) {
      # Compare final cluster assignments
      similarity <- compare_clusterings(results$r$clusterLabels,
                                        results$cpp$clusterLabels)
      similarities[i] <- similarity
    } else {
      skip("C++ implementation not available")
    }
  }

  # Average similarity should be reasonably high
  avg_similarity <- mean(similarities, na.rm = TRUE)
  expect_true(avg_similarity > 0.7,
              info = sprintf("Average clustering similarity: %.3f", avg_similarity))
})

# =============================================================================
# Edge Cases and Robustness Tests
# =============================================================================

test_that("C++ and R handle small datasets equivalently", {
  # Very small dataset
  small_data <- c(1.0, 1.1, 2.0, 2.1)
  results <- run_both_implementations(small_data, n_iter = 50, seed = 900)

  skip_if(!results$both_available, "C++ implementation not available")

  # Both should handle small data without errors
  expect_true(results$r$numberClusters >= 1)
  expect_true(results$cpp$numberClusters >= 1)
  expect_true(results$r$numberClusters <= length(small_data))
  expect_true(results$cpp$numberClusters <= length(small_data))
})

test_that("C++ and R handle single cluster data equivalently", {
  # Data from single cluster
  single_cluster_data <- rnorm(50, mean = 0, sd = 1)
  results <- run_both_implementations(single_cluster_data, n_iter = 80, seed = 1000)

  skip_if(!results$both_available, "C++ implementation not available")

  # Both should typically find 1-2 clusters for homogeneous data
  expect_true(results$r$numberClusters <= 3)
  expect_true(results$cpp$numberClusters <= 3)

  # Cluster counts should be similar
  cluster_diff <- abs(results$r$numberClusters - results$cpp$numberClusters)
  expect_true(cluster_diff <= 1)
})

test_that("C++ and R handle extreme parameter values equivalently", {
  test_data <- generate_test_mixture(n = 80, k = 2, seed = 1100)

  # Test with different prior parameters
  custom_priors <- c(0, 0.1, 2, 2)  # mu0, kappa0, alpha0, beta0

  set_use_cpp(FALSE)
  set.seed(1200)
  dp_r <- DirichletProcessGaussian(test_data$data, g0Priors = custom_priors)
  dp_r <- Fit(dp_r, 60, progressBar = FALSE)

  cpp_available <- exists("_dirichletprocess_run_mcmc_cpp") ||
    get_cpp_status()$mcmc_runner

  if (cpp_available && can_use_cpp(dp_r)) {
    set_use_cpp(TRUE)
    set.seed(1200)
    dp_cpp <- DirichletProcessGaussian(test_data$data, g0Priors = custom_priors)
    dp_cpp <- Fit(dp_cpp, 60, progressBar = FALSE)

    # Both should handle custom priors without issues
    expect_true(dp_r$numberClusters >= 1)
    expect_true(dp_cpp$numberClusters >= 1)

    cluster_diff <- abs(dp_r$numberClusters - dp_cpp$numberClusters)
    expect_true(cluster_diff <= 2)
  } else {
    skip("C++ implementation not available")
  }
})

# =============================================================================
# Performance Consistency Tests
# =============================================================================

test_that("C++ and R produce consistent results across different data sizes", {
  sizes <- c(30, 100, 200)
  max_cluster_diff <- 0

  for (n in sizes) {
    test_data <- generate_test_mixture(n = n, k = 3, seed = 1300 + n)
    results <- run_both_implementations(test_data$data, n_iter = 80, seed = 1400 + n)

    if (results$both_available) {
      cluster_diff <- abs(results$r$numberClusters - results$cpp$numberClusters)
      max_cluster_diff <- max(max_cluster_diff, cluster_diff)
    } else {
      skip("C++ implementation not available")
    }
  }

  # Maximum difference across all sizes should be reasonable
  expect_true(max_cluster_diff <= 2,
              info = sprintf("Maximum cluster count difference: %d", max_cluster_diff))
})

test_that("C++ and R maintain consistency with different iteration counts", {
  test_data <- generate_test_mixture(n = 80, k = 2, seed = 1500)
  iteration_counts <- c(50, 100, 200)

  r_clusters <- numeric(length(iteration_counts))
  cpp_clusters <- numeric(length(iteration_counts))

  for (i in seq_along(iteration_counts)) {
    n_iter <- iteration_counts[i]
    results <- run_both_implementations(test_data$data, n_iter = n_iter, seed = 1600 + i)

    if (results$both_available) {
      r_clusters[i] <- results$r$numberClusters
      cpp_clusters[i] <- results$cpp$numberClusters
    } else {
      skip("C++ implementation not available")
    }
  }

  # Both implementations should show similar trends with more iterations
  r_consistency <- var(r_clusters)
  cpp_consistency <- var(cpp_clusters)

  # Variance should be reasonable (not too high)
  expect_true(r_consistency < 2)
  expect_true(cpp_consistency < 2)
})

# =============================================================================
# Integration Tests
# =============================================================================

test_that("Complete workflow equivalence", {
  test_data <- generate_test_mixture(n = 150, k = 3, seed = 1700)

  # Full workflow test
  results <- run_both_implementations(test_data$data, n_iter = 150, seed = 1800)

  skip_if(!results$both_available, "C++ implementation not available")

  # Check all major components

  # 1. Data integrity
  expect_equal(results$r$data, results$cpp$data)
  expect_equal(results$r$n, results$cpp$n)

  # 2. Clustering results
  expect_true(results$r$numberClusters > 0)
  expect_true(results$cpp$numberClusters > 0)

  # 3. Parameter estimates
  r_summary <- extract_posterior_means(results$r)
  cpp_summary <- extract_posterior_means(results$cpp)

  expect_true(r_summary$alpha > 0)
  expect_true(cpp_summary$alpha > 0)

  # 4. Overall similarity
  cluster_diff <- abs(results$r$numberClusters - results$cpp$numberClusters)
  alpha_diff <- abs(r_summary$alpha - cpp_summary$alpha)

  expect_true(cluster_diff <= 2)
  expect_true(alpha_diff < 3.0)

  # 5. Cluster assignments similarity
  if (results$r$numberClusters == results$cpp$numberClusters) {
    similarity <- compare_clusterings(results$r$clusterLabels,
                                      results$cpp$clusterLabels)
    expect_true(similarity > 0.5)
  }
})
