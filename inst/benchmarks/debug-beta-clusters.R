# Debug script for Beta Distribution test failures
# This script will help identify the root cause of pointsPerCluster mismatches

library(dirichletprocess)
library(testthat)

# Enable detailed output
options(dirichletprocess.debug = TRUE)

# Create a debug wrapper function to track state changes
debug_beta_dp <- function(data, n_iter = 100, verbose = TRUE) {

  cat("\n========== BETA DP DEBUG SESSION ==========\n")
  cat("Data points:", length(data), "\n")
  cat("Data range: [", min(data), ",", max(data), "]\n\n")

  # Initialize the DP
  dp <- DirichletProcessBeta(data)

  # Check initial state
  cat("=== INITIAL STATE ===\n")
  cat("Number of clusters:", dp$numberClusters, "\n")
  cat("Cluster labels:", paste(dp$clusterLabels, collapse = ", "), "\n")
  cat("Points per cluster:", paste(dp$pointsPerCluster, collapse = ", "), "\n")
  cat("Sum of pointsPerCluster:", sum(dp$pointsPerCluster), "\n")
  cat("Expected sum (n):", length(data), "\n")
  cat("Initial mismatch:", sum(dp$pointsPerCluster) - length(data), "\n\n")

  # Track state through iterations
  state_history <- list()

  # Custom update function with debugging
  debug_update <- function(dp_obj, iteration) {

    # Store pre-update state
    pre_state <- list(
      iteration = iteration,
      n_clusters = dp_obj$numberClusters,
      cluster_labels = dp_obj$clusterLabels,
      points_per_cluster = dp_obj$pointsPerCluster,
      sum_points = sum(dp_obj$pointsPerCluster),
      alpha = dp_obj$alpha
    )

    # Perform cluster component update
    tryCatch({
      # Update cluster assignments
      for (i in seq_along(dp_obj$data)) {
        old_cluster <- dp_obj$clusterLabels[i]

        # Remove point from old cluster
        if (old_cluster > 0 && old_cluster <= length(dp_obj$pointsPerCluster)) {
          dp_obj$pointsPerCluster[old_cluster] <- dp_obj$pointsPerCluster[old_cluster] - 1
        }

        # Calculate probabilities for each cluster
        cluster_probs <- numeric(dp_obj$numberClusters + 1)

        # Existing clusters
        for (j in 1:dp_obj$numberClusters) {
          if (dp_obj$pointsPerCluster[j] > 0 || j == old_cluster) {
            n_j <- dp_obj$pointsPerCluster[j]
            if (j == old_cluster && n_j == 0) n_j <- 0  # Handle empty cluster

            # Calculate likelihood
            lik <- Likelihood(dp_obj$mixingDistribution,
                              matrix(dp_obj$data[i], ncol = 1),
                              dp_obj$clusterParameters[[j]])

            cluster_probs[j] <- n_j * lik
          }
        }

        # New cluster probability
        new_cluster_param <- PriorDraw(dp_obj$mixingDistribution)
        new_lik <- Likelihood(dp_obj$mixingDistribution,
                              matrix(dp_obj$data[i], ncol = 1),
                              new_cluster_param)
        cluster_probs[dp_obj$numberClusters + 1] <- dp_obj$alpha * new_lik

        # Normalize probabilities
        cluster_probs <- cluster_probs / sum(cluster_probs)

        # Sample new cluster
        new_cluster <- sample.int(dp_obj$numberClusters + 1, 1, prob = cluster_probs)

        # Assign to cluster
        if (new_cluster <= dp_obj$numberClusters) {
          dp_obj$clusterLabels[i] <- new_cluster
          dp_obj$pointsPerCluster[new_cluster] <- dp_obj$pointsPerCluster[new_cluster] + 1
        } else {
          # Create new cluster
          dp_obj$numberClusters <- dp_obj$numberClusters + 1
          dp_obj$clusterLabels[i] <- dp_obj$numberClusters
          dp_obj$pointsPerCluster <- c(dp_obj$pointsPerCluster, 1)
          dp_obj$clusterParameters[[dp_obj$numberClusters]] <- new_cluster_param
        }
      }

      # Clean up empty clusters
      non_empty <- which(dp_obj$pointsPerCluster > 0)
      if (length(non_empty) < dp_obj$numberClusters) {
        # Remap cluster labels
        new_labels <- integer(length(dp_obj$clusterLabels))
        new_params <- list()

        for (j in seq_along(non_empty)) {
          old_idx <- non_empty[j]
          new_labels[dp_obj$clusterLabels == old_idx] <- j
          new_params[[j]] <- dp_obj$clusterParameters[[old_idx]]
        }

        dp_obj$clusterLabels <- new_labels
        dp_obj$clusterParameters <- new_params
        dp_obj$pointsPerCluster <- dp_obj$pointsPerCluster[non_empty]
        dp_obj$numberClusters <- length(non_empty)
      }

    }, error = function(e) {
      cat("\nERROR in cluster update at iteration", iteration, ":\n")
      cat(conditionMessage(e), "\n")
      print(traceback())
    })

    # Store post-update state
    post_state <- list(
      iteration = iteration,
      n_clusters = dp_obj$numberClusters,
      cluster_labels = dp_obj$clusterLabels,
      points_per_cluster = dp_obj$pointsPerCluster,
      sum_points = sum(dp_obj$pointsPerCluster),
      alpha = dp_obj$alpha
    )

    # Check for issues
    issues <- character()

    if (sum(dp_obj$pointsPerCluster) != length(dp_obj$data)) {
      issues <- c(issues, sprintf("Points mismatch: sum=%d, n=%d",
                                  sum(dp_obj$pointsPerCluster),
                                  length(dp_obj$data)))
    }

    if (any(dp_obj$clusterLabels <= 0)) {
      issues <- c(issues, "Found non-positive cluster labels")
    }

    if (any(dp_obj$clusterLabels > dp_obj$numberClusters)) {
      issues <- c(issues, "Found cluster labels exceeding numberClusters")
    }

    if (length(dp_obj$pointsPerCluster) != dp_obj$numberClusters) {
      issues <- c(issues, sprintf("Length mismatch: pointsPerCluster=%d, numberClusters=%d",
                                  length(dp_obj$pointsPerCluster),
                                  dp_obj$numberClusters))
    }

    # Report issues
    if (length(issues) > 0 && verbose) {
      cat("\n!!! ISSUES DETECTED at iteration", iteration, "!!!\n")
      for (issue in issues) {
        cat("  -", issue, "\n")
      }
      cat("\nState details:\n")
      cat("  Cluster labels:", paste(dp_obj$clusterLabels, collapse = ", "), "\n")
      cat("  Points per cluster:", paste(dp_obj$pointsPerCluster, collapse = ", "), "\n")
      cat("  Unique clusters in labels:", paste(sort(unique(dp_obj$clusterLabels)), collapse = ", "), "\n")
    }

    # Store state change
    state_history[[length(state_history) + 1]] <- list(
      pre = pre_state,
      post = post_state,
      issues = issues
    )

    # Update cluster parameters
    dp_obj <- ClusterParameterUpdate(dp_obj)

    # Update alpha
    dp_obj <- UpdateAlpha(dp_obj)

    return(dp_obj)
  }

  # Run iterations with debugging
  cat("\n=== RUNNING ITERATIONS ===\n")
  for (iter in 1:n_iter) {
    if (iter %% 10 == 0 && verbose) {
      cat("Iteration", iter, "- Clusters:", dp$numberClusters,
          "- Sum points:", sum(dp$pointsPerCluster), "\n")
    }

    dp <- debug_update(dp, iter)

    # Stop if critical error
    if (sum(dp$pointsPerCluster) != length(dp$data)) {
      cat("\n!!! STOPPING: Critical mismatch detected !!!\n")
      break
    }
  }

  # Final state check
  cat("\n=== FINAL STATE ===\n")
  cat("Number of clusters:", dp$numberClusters, "\n")
  cat("Cluster labels:", paste(head(dp$clusterLabels, 20), collapse = ", "),
      if(length(dp$clusterLabels) > 20) "...", "\n")
  cat("Points per cluster:", paste(dp$pointsPerCluster, collapse = ", "), "\n")
  cat("Sum of pointsPerCluster:", sum(dp$pointsPerCluster), "\n")
  cat("Expected sum (n):", length(dp$data), "\n")
  cat("Final mismatch:", sum(dp$pointsPerCluster) - length(dp$data), "\n")

  # Analyze state history for patterns
  if (length(state_history) > 0) {
    cat("\n=== ISSUE SUMMARY ===\n")
    all_issues <- unlist(lapply(state_history, function(x) x$issues))
    if (length(all_issues) > 0) {
      issue_table <- table(all_issues)
      for (i in 1:length(issue_table)) {
        cat(names(issue_table)[i], ":", issue_table[i], "occurrences\n")
      }
    } else {
      cat("No issues detected during iterations.\n")
    }
  }

  return(list(
    dp = dp,
    state_history = state_history,
    initial_n = length(data),
    final_sum = sum(dp$pointsPerCluster),
    mismatch = sum(dp$pointsPerCluster) - length(data)
  ))
}

# Test with different scenarios
cat("\n############################################\n")
cat("# TEST 1: Small dataset (n=10)\n")
cat("############################################\n")
set.seed(123)
result1 <- debug_beta_dp(rbeta(10, 2, 5), n_iter = 20, verbose = TRUE)

cat("\n############################################\n")
cat("# TEST 2: Two-cluster data (n=20)\n")
cat("############################################\n")
set.seed(456)
data2 <- c(rbeta(10, 2, 8), rbeta(10, 8, 2))
result2 <- debug_beta_dp(data2, n_iter = 20, verbose = TRUE)

cat("\n############################################\n")
cat("# TEST 3: Larger dataset (n=50)\n")
cat("############################################\n")
set.seed(789)
result3 <- debug_beta_dp(rbeta(50, 3, 3), n_iter = 20, verbose = TRUE)

# Additional diagnostic function
diagnose_cluster_state <- function(dp) {
  cat("\n=== DETAILED CLUSTER DIAGNOSTICS ===\n")

  # Check each data point
  cat("\nData point analysis:\n")
  label_counts <- table(dp$clusterLabels)
  cat("Cluster label counts:\n")
  print(label_counts)

  cat("\nExpected vs Actual points per cluster:\n")
  for (i in 1:dp$numberClusters) {
    expected <- sum(dp$clusterLabels == i)
    actual <- dp$pointsPerCluster[i]
    if (expected != actual) {
      cat(sprintf("  Cluster %d: expected=%d, actual=%d [MISMATCH]\n",
                  i, expected, actual))
    } else {
      cat(sprintf("  Cluster %d: %d points [OK]\n", i, actual))
    }
  }

  # Check for orphaned labels
  orphaned <- setdiff(unique(dp$clusterLabels), 1:dp$numberClusters)
  if (length(orphaned) > 0) {
    cat("\nORPHANED LABELS FOUND:", paste(orphaned, collapse = ", "), "\n")
  }

  # Check for missing labels
  missing <- setdiff(1:dp$numberClusters, unique(dp$clusterLabels))
  if (length(missing) > 0) {
    cat("\nMISSING LABELS (empty clusters):", paste(missing, collapse = ", "), "\n")
  }
}

# Run diagnostics on failed cases
if (result3$mismatch != 0) {
  diagnose_cluster_state(result3$dp)
}

# Test the specific failing test case
cat("\n############################################\n")
cat("# TEST 4: Reproduce test failure scenario\n")
cat("############################################\n")
set.seed(42)  # Use same seed as test
dp_test <- DirichletProcessBeta(rbeta(20, 2, 5))
dp_test <- Initialise(dp_test)

# Run a few updates manually
for (i in 1:5) {
  cat("\n--- Update", i, "---\n")
  cat("Before: clusters =", dp_test$numberClusters,
      ", sum =", sum(dp_test$pointsPerCluster), "\n")

  dp_test <- ClusterComponentUpdate(dp_test)

  cat("After: clusters =", dp_test$numberClusters,
      ", sum =", sum(dp_test$pointsPerCluster), "\n")

  if (sum(dp_test$pointsPerCluster) != length(dp_test$data)) {
    cat("MISMATCH DETECTED!\n")
    diagnose_cluster_state(dp_test)
    break
  }
}
