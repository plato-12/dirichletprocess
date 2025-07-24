# tests/integration/stress_tests.R

run_stress_tests <- function() {
  stress_results <- list()

  # 1. Large dataset test
  cat("\n1. Testing large dataset (n=50,000)...\n")
  large_data <- rnorm(50000)

  tryCatch({
    dp <- DirichletProcessGaussian(large_data)
    dp <- Fit(dp, its = 10)
    stress_results$large_dataset <- "PASSED"
  }, error = function(e) {
    stress_results$large_dataset <- paste("FAILED:", e$message)
  })

  # 2. Many clusters test
  cat("\n2. Testing many clusters scenario...\n")
  many_clusters_data <- c()
  for (i in 1:50) {
    many_clusters_data <- c(many_clusters_data, rnorm(20, mean = i * 5))
  }

  tryCatch({
    dp <- DirichletProcessGaussian(many_clusters_data)
    dp <- Fit(dp, its = 100)
    stress_results$many_clusters <- paste("PASSED - Found", dp$numberClusters, "clusters")
  }, error = function(e) {
    stress_results$many_clusters <- paste("FAILED:", e$message)
  })

  # 3. High dimension test (MVNormal)
  cat("\n3. Testing high dimensions (d=20)...\n")
  high_dim_data <- matrix(rnorm(1000 * 20), ncol = 20)

  tryCatch({
    dp <- DirichletProcessMvnormal(high_dim_data)
    dp <- Fit(dp, its = 10)
    stress_results$high_dimension <- "PASSED"
  }, error = function(e) {
    stress_results$high_dimension <- paste("FAILED:", e$message)
  })

  # 4. Extreme parameter values
  cat("\n4. Testing extreme parameters...\n")
  extreme_data <- c(rnorm(100, sd = 0.01), rnorm(100, sd = 1000))

  tryCatch({
    dp <- DirichletProcessGaussian(extreme_data)
    dp <- Fit(dp, its = 50)
    stress_results$extreme_params <- "PASSED"
  }, error = function(e) {
    stress_results$extreme_params <- paste("FAILED:", e$message)
  })

  # 5. Concurrent execution test
  cat("\n5. Testing concurrent execution...\n")
  if (requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(4)

    tryCatch({
      parallel::clusterEvalQ(cl, {
        library(dirichletprocess)
        set_use_cpp(TRUE)
      })

      results <- parallel::parLapply(cl, 1:4, function(i) {
        dp <- DirichletProcessGaussian(rnorm(1000))
        dp <- Fit(dp, its = 100)
        return(dp$numberClusters)
      })

      stress_results$concurrent <- "PASSED"
    }, error = function(e) {
      stress_results$concurrent <- paste("FAILED:", e$message)
    }, finally = {
      parallel::stopCluster(cl)
    })
  }

  return(stress_results)
}

# Extended stress tests
run_extended_stress_tests <- function() {
  cat("\n=== EXTENDED STRESS TESTS ===\n")

  results <- list()

  # 1. Imbalanced clusters
  cat("\n1. Testing highly imbalanced clusters...\n")
  imbalanced_data <- c(rnorm(1000, mean = 0),
                       rnorm(10, mean = 10),
                       rnorm(5, mean = -10))

  results$imbalanced <- tryCatch({
    dp <- DirichletProcessGaussian(imbalanced_data)
    dp <- Fit(dp, its = 100)
    list(
      status = "PASSED",
      clusters_found = dp$numberClusters,
      cluster_sizes = table(dp$clusterLabels)
    )
  }, error = function(e) {
    list(status = "FAILED", error = e$message)
  })

  # 2. Numerical precision stress
  cat("\n2. Testing numerical precision limits...\n")

  # Very small scale
  small_scale_data <- rnorm(100) * 1e-15
  results$small_scale <- tryCatch({
    dp <- DirichletProcessGaussian(small_scale_data)
    dp <- Fit(dp, its = 50)
    "PASSED"
  }, error = function(e) {
    paste("FAILED:", e$message)
  })

  # Very large scale
  large_scale_data <- rnorm(100) * 1e15
  results$large_scale <- tryCatch({
    dp <- DirichletProcessGaussian(large_scale_data)
    dp <- Fit(dp, its = 50)
    "PASSED"
  }, error = function(e) {
    paste("FAILED:", e$message)
  })

  # 3. Long running test
  cat("\n3. Testing long running MCMC (10,000 iterations)...\n")
  long_run_data <- generate_test_data("normal", 500)

  results$long_run <- tryCatch({
    start_time <- Sys.time()
    dp <- DirichletProcessGaussian(long_run_data)
    dp <- Fit(dp, its = 10000)
    end_time <- Sys.time()

    list(
      status = "PASSED",
      runtime = as.numeric(end_time - start_time, units = "secs"),
      final_clusters = dp$numberClusters,
      alpha_mean = mean(dp$alphaChain)
    )
  }, error = function(e) {
    list(status = "FAILED", error = e$message)
  })

  # 4. Rapid cluster changes
  cat("\n4. Testing rapid cluster changes...\n")
  # Data that encourages cluster splitting/merging
  dynamic_data <- c()
  for (i in 1:10) {
    dynamic_data <- c(dynamic_data, rnorm(50, mean = sin(i) * 10))
  }

  results$dynamic_clusters <- tryCatch({
    dp <- DirichletProcessGaussian(dynamic_data)
    dp <- Fit(dp, its = 500)

    # Analyze cluster stability
    cluster_counts <- sapply(dp$labelsChain, function(x) length(unique(x)))

    list(
      status = "PASSED",
      min_clusters = min(cluster_counts),
      max_clusters = max(cluster_counts),
      cluster_variance = var(cluster_counts)
    )
  }, error = function(e) {
    list(status = "FAILED", error = e$message)
  })

  # 5. Mixed distribution types stress test
  cat("\n5. Testing all distributions with edge cases...\n")

  edge_case_data <- list(
    normal = c(rnorm(50, -100), rnorm(50, 100)),
    exponential = c(rexp(50, 0.01), rexp(50, 100)),
    beta = c(rbeta(50, 0.1, 0.1), rbeta(50, 10, 10)),
    weibull = c(rweibull(50, 0.1), rweibull(50, 10)),
    mvnormal = rbind(
      mvtnorm::rmvnorm(50, rep(-10, 5), diag(5)),
      mvtnorm::rmvnorm(50, rep(10, 5), diag(5))
    )
  )

  for (dist in names(edge_case_data)) {
    cat("    Testing", dist, "...")

    results[[paste0(dist, "_edge")]] <- tryCatch({
      dp <- create_dp_object(dist, edge_case_data[[dist]])
      dp <- Fit(dp, its = 100)
      "PASSED"
    }, error = function(e) {
      paste("FAILED:", e$message)
    })

    cat(" ", results[[paste0(dist, "_edge")]], "\n")
  }

  return(results)
}

# Robustness tests
test_robustness <- function() {
  cat("\n=== ROBUSTNESS TESTS ===\n")

  results <- list()

  # 1. Repeated fitting stability
  cat("\n1. Testing repeated fitting stability...\n")
  test_data <- generate_test_data("normal", 200)

  cluster_results <- numeric(10)
  alpha_results <- numeric(10)

  for (i in 1:10) {
    set.seed(123)  # Same seed each time
    dp <- DirichletProcessGaussian(test_data)
    dp <- Fit(dp, its = 100)

    cluster_results[i] <- dp$numberClusters
    alpha_results[i] <- dp$alpha
  }

  results$repeated_fitting <- list(
    cluster_variance = var(cluster_results),
    alpha_variance = var(alpha_results),
    all_identical = length(unique(cluster_results)) == 1
  )

  # 2. Recovery from bad initialization
  cat("\n2. Testing recovery from bad initialization...\n")
  test_data <- generate_test_data("exponential", 300)
  dp <- DirichletProcessExponential(test_data)

  # Force bad initialization
  dp$clusterLabels <- rep(1, length(test_data))  # All in one cluster
  dp$numberClusters <- 1
  dp$alpha <- 0.001  # Very low alpha

  results$bad_init_recovery <- tryCatch({
    dp <- Fit(dp, its = 200)
    list(
      status = "PASSED",
      final_clusters = dp$numberClusters,
      recovered = dp$numberClusters > 1
    )
  }, error = function(e) {
    list(status = "FAILED", error = e$message)
  })

  # 3. Handling of outliers
  cat("\n3. Testing outlier handling...\n")
  # Normal data with outliers
  outlier_data <- c(rnorm(200),
                    runif(10, min = -50, max = -40),  # Left outliers
                    runif(10, min = 40, max = 50))     # Right outliers

  results$outlier_handling <- tryCatch({
    dp <- DirichletProcessGaussian(outlier_data)
    dp <- Fit(dp, its = 200)

    # Check if outliers are in separate clusters
    outlier_indices <- c(201:220)
    outlier_clusters <- dp$clusterLabels[outlier_indices]
    main_clusters <- dp$clusterLabels[1:200]

    list(
      status = "PASSED",
      total_clusters = dp$numberClusters,
      outliers_separated = length(intersect(unique(outlier_clusters),
                                            unique(main_clusters))) == 0
    )
  }, error = function(e) {
    list(status = "FAILED", error = e$message)
  })

  return(results)
}

# Performance under extreme conditions
test_extreme_performance <- function() {
  cat("\n=== EXTREME PERFORMANCE TESTS ===\n")

  results <- list()

  # 1. Maximum practical dataset size
  cat("\n1. Finding maximum practical dataset size...\n")

  sizes <- c(10000, 25000, 50000, 75000, 100000)
  max_size <- 0

  for (size in sizes) {
    cat("  Trying n =", size, "...")

    result <- tryCatch({
      test_data <- rnorm(size)
      start_time <- Sys.time()

      dp <- DirichletProcessGaussian(test_data)
      dp <- Fit(dp, its = 5)  # Just 5 iterations

      end_time <- Sys.time()
      runtime <- as.numeric(end_time - start_time, units = "secs")

      if (runtime < 60) {  # Less than 1 minute
        max_size <- size
        cat(" SUCCESS (", round(runtime, 1), "s)\n")
        TRUE
      } else {
        cat(" Too slow (", round(runtime, 1), "s)\n")
        FALSE
      }
    }, error = function(e) {
      cat(" FAILED\n")
      FALSE
    })

    if (!result) break
  }

  results$max_dataset_size <- max_size

  # 2. Maximum dimensions for multivariate
  cat("\n2. Finding maximum dimensions for MVNormal...\n")

  dimensions <- c(10, 25, 50, 75, 100)
  max_dim <- 0

  for (d in dimensions) {
    cat("  Trying d =", d, "...")

    result <- tryCatch({
      test_data <- matrix(rnorm(1000 * d), ncol = d)

      dp <- DirichletProcessMvnormal(test_data)
      dp <- Fit(dp, its = 10)

      max_dim <- d
      cat(" SUCCESS\n")
      TRUE
    }, error = function(e) {
      cat(" FAILED\n")
      FALSE
    })

    if (!result) break
  }

  results$max_dimensions <- max_dim

  return(results)
}

# Run comprehensive stress test suite
run_comprehensive_stress_tests <- function() {
  all_results <- list()

  cat("\n========================================\n")
  cat("  COMPREHENSIVE STRESS TEST SUITE\n")
  cat("========================================\n")

  # Basic stress tests
  all_results$basic <- run_stress_tests()

  # Extended stress tests
  all_results$extended <- run_extended_stress_tests()

  # Robustness tests
  all_results$robustness <- test_robustness()

  # Extreme performance tests
  all_results$extreme <- test_extreme_performance()

  # Generate stress test report
  generate_stress_report(all_results)

  return(all_results)
}

generate_stress_report <- function(results) {
  cat("\n\n=== STRESS TEST REPORT ===\n")
  cat("Generated:", format(Sys.time()), "\n\n")

  # Basic stress tests
  cat("Basic Stress Tests:\n")
  for (test in names(results$basic)) {
    cat("  ", test, ":", results$basic[[test]], "\n")
  }

  # Extended tests summary
  cat("\nExtended Tests:\n")
  if (!is.null(results$extended$long_run$runtime)) {
    cat("  Long run (10k iterations):",
        round(results$extended$long_run$runtime, 1), "seconds\n")
  }

  # Robustness summary
  cat("\nRobustness:\n")
  if (!is.null(results$robustness$repeated_fitting)) {
    cat("  Repeated fitting stable:",
        results$robustness$repeated_fitting$all_identical, "\n")
  }

  # Extreme performance
  cat("\nExtreme Performance:\n")
  if (!is.null(results$extreme$max_dataset_size)) {
    cat("  Max dataset size:", results$extreme$max_dataset_size, "\n")
  }
  if (!is.null(results$extreme$max_dimensions)) {
    cat("  Max dimensions:", results$extreme$max_dimensions, "\n")
  }

  # Save detailed results
  saveRDS(results, "stress_test_results.rds")
  cat("\nDetailed results saved to: stress_test_results.rds\n")
}
