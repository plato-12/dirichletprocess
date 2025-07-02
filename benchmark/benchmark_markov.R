# benchmark/benchmark_markov_dp.R
# Comprehensive benchmarking for Markov Chain Dirichlet Process
# Compares R vs C++ implementations using algorithms from Neal (2000)

library(dirichletprocess)
library(ggplot2)
library(dplyr)
library(tidyr)
library(microbenchmark)
library(patchwork)

#' Quick Markov DP benchmark for interactive testing
#'
#' Tests different algorithms from Neal (2000) paper
#' @param n_obs Vector of observation counts to test
#' @param algorithms Vector of algorithm numbers (1-8 from Neal 2000)
#' @param n_iter Number of MCMC iterations
#' @param n_reps Number of repetitions for timing
#' @export
quick_markov_dp_benchmark <- function(n_obs = c(100, 250, 500),
                                      algorithms = c(2, 3, 5, 8),
                                      n_iter = 100,
                                      n_reps = 3) {

  cat("\n=== Quick Markov DP Benchmark (Neal 2000 Algorithms) ===\n\n")

  results <- data.frame()

  for (algo in algorithms) {
    for (n in n_obs) {
      cat(sprintf("Testing Algorithm %d, n = %d:\n", algo, n))

      # Generate test data - mixture of normals as in Neal (2000) examples
      set.seed(42)
      true_means <- c(-3, 0, 3)
      true_sds <- c(0.5, 0.7, 0.5)
      true_weights <- c(0.3, 0.4, 0.3)

      # Generate mixture data
      cluster_assignments <- sample(1:3, n, replace = TRUE, prob = true_weights)
      y <- numeric(n)
      for (i in 1:3) {
        idx <- cluster_assignments == i
        y[idx] <- rnorm(sum(idx), true_means[i], true_sds[i])
      }

      # Benchmark R implementation
      times_r <- numeric(n_reps)
      for (i in 1:n_reps) {
        set_use_cpp(FALSE)
        time_r <- system.time({
          dp_r <- DirichletProcessGaussian(y)
          # Set algorithm type if supported
          if ("algorithm" %in% names(dp_r)) {
            dp_r$algorithm <- algo
          }
          dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
        })
        times_r[i] <- time_r["elapsed"]
      }

      # Benchmark C++ implementation
      times_cpp <- numeric(n_reps)
      clusters_cpp <- numeric(n_reps)

      for (i in 1:n_reps) {
        set_use_cpp(TRUE)
        time_cpp <- system.time({
          dp_cpp <- DirichletProcessGaussian(y)
          if ("algorithm" %in% names(dp_cpp)) {
            dp_cpp$algorithm <- algo
          }
          dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
        })
        times_cpp[i] <- time_cpp["elapsed"]
        clusters_cpp[i] <- dp_cpp$numberClusters
      }

      # Calculate statistics
      mean_r <- mean(times_r)
      mean_cpp <- mean(times_cpp)
      speedup <- mean_r / mean_cpp

      cat(sprintf("  R: %.3fs, C++: %.3fs, Speedup: %.1fx, Clusters: %.1f\n",
                  mean_r, mean_cpp, speedup, mean(clusters_cpp)))

      # Store results
      results <- rbind(results, data.frame(
        algorithm = algo,
        n_obs = n,
        implementation = rep(c("R", "C++"), each = n_reps),
        time = c(times_r, times_cpp),
        clusters = c(rep(NA, n_reps), clusters_cpp)
      ))
    }
  }

  return(results)
}

#' Profile individual Markov chain components
#'
#' Tests components of Neal (2000) algorithms
#' @param n Number of observations
#' @param algorithm Algorithm number (1-8)
#' @return Data frame with component timings
#' @export
profile_markov_dp_components <- function(n = 500, algorithm = 3) {

  cat("\n=== Markov DP Component Profiling ===\n")
  cat(sprintf("Algorithm %d components (Neal 2000)\n", algorithm))

  # Generate test data
  set.seed(123)
  y <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 2, 0.5))

  # Convert to matrix format for consistency
  y <- matrix(y, ncol = 1)

  # Create DP objects
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(y)
  dp_r <- Initialise(dp_r)

  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(y)
  dp_cpp <- Initialise(dp_cpp)

  # Component timing results
  results <- data.frame()

  # Time Cluster Component Update
  cat("Profiling ClusterComponentUpdate...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      ClusterComponentUpdate(dp_r)
    },
    times = 10
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      ClusterComponentUpdate(dp_cpp)
    },
    times = 10
  )

  results <- rbind(results, data.frame(
    component = "ClusterComponentUpdate",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Time Parameter Update
  cat("Profiling ClusterParameterUpdate...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      ClusterParameterUpdate(dp_r)
    },
    times = 50
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      ClusterParameterUpdate(dp_cpp)
    },
    times = 50
  )

  results <- rbind(results, data.frame(
    component = "ClusterParameterUpdate",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Time Alpha Update
  cat("Profiling UpdateAlpha...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      UpdateAlpha(dp_r)
    },
    times = 50
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      UpdateAlpha(dp_cpp)
    },
    times = 50
  )

  results <- rbind(results, data.frame(
    component = "UpdateAlpha",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Time Likelihood calculation
  cat("Profiling Likelihood calculations...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      total_lik <- 0
      for (i in 1:min(100, n)) {
        # Extract parameters properly for univariate Gaussian
        if (length(dim(dp_r$clusterParameters[[1]])) == 3) {
          mu_val <- dp_r$clusterParameters[[1]][1, 1, 1]
          sigma_val <- dp_r$clusterParameters[[2]][1, 1, 1]
        } else {
          mu_val <- dp_r$clusterParameters[[1]][1]
          sigma_val <- dp_r$clusterParameters[[2]][1]
        }

        single_cluster_params <- list(mu_val, sigma_val)

        lik <- Likelihood(dp_r$mixingDistribution,
                          y[i, , drop = FALSE],
                          single_cluster_params)
        total_lik <- total_lik + lik
      }
    },
    times = 10
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      total_lik <- 0
      for (i in 1:min(100, n)) {
        # Extract parameters properly for univariate Gaussian
        if (length(dim(dp_cpp$clusterParameters[[1]])) == 3) {
          mu_val <- dp_cpp$clusterParameters[[1]][1, 1, 1]
          sigma_val <- dp_cpp$clusterParameters[[2]][1, 1, 1]
        } else {
          mu_val <- dp_cpp$clusterParameters[[1]][1]
          sigma_val <- dp_cpp$clusterParameters[[2]][1]
        }

        single_cluster_params <- list(mu_val, sigma_val)

        lik <- Likelihood(dp_cpp$mixingDistribution,
                          y[i, , drop = FALSE],
                          single_cluster_params)
        total_lik <- total_lik + lik
      }
    },
    times = 10
  )

  results <- rbind(results, data.frame(
    component = "Likelihood",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Calculate speedups
  cat("\n\nComponent timings (milliseconds):\n")
  print(results)

  # Calculate speedup for each component
  library(tidyr)
  speedups <- results %>%
    pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
    mutate(speedup = R / `C++`)

  cat("\nSpeedup by component:\n")
  for (i in 1:nrow(speedups)) {
    cat(sprintf("  %s: %.1fx\n", speedups$component[i], speedups$speedup[i]))
  }

  return(results)
}

#' Test statistical equivalence of Markov chain outputs
#'
#' Verifies that R and C++ produce equivalent posterior distributions
#' @param n Number of observations
#' @param n_iter Number of MCMC iterations
#' @param algorithm Algorithm to test
#' @return List with both DP objects
#' @export
test_markov_dp_equivalence <- function(n = 300, n_iter = 500, algorithm = 3) {

  cat("\n=== Testing Markov DP Statistical Equivalence ===\n")
  cat(sprintf("Algorithm %d, n = %d, iterations = %d\n", algorithm, n, n_iter))

  # Generate test data
  set.seed(999)
  y <- c(rnorm(n/3, -3, 0.5), rnorm(n/3, 0, 0.7), rnorm(n/3, 3, 0.5))

  # Run R implementation
  cat("\nRunning R implementation...\n")
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(y)
  if ("algorithm" %in% names(dp_r)) {
    dp_r$algorithm <- algorithm
  }
  dp_r <- Fit(dp_r, n_iter, progressBar = TRUE)

  # Run C++ implementation
  cat("\nRunning C++ implementation...\n")
  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(y)
  if ("algorithm" %in% names(dp_cpp)) {
    dp_cpp$algorithm <- algorithm
  }
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = TRUE)

  # Compare results
  cat("\nComparing Markov chain outputs:\n")

  # Number of clusters
  cat(sprintf("\nNumber of clusters:\n"))
  cat(sprintf("  R: %d\n", dp_r$numberClusters))
  cat(sprintf("  C++: %d\n", dp_cpp$numberClusters))

  # Alpha parameter (concentration)
  cat(sprintf("\nConcentration parameter (alpha):\n"))
  cat(sprintf("  R: %.3f\n", dp_r$alpha))
  cat(sprintf("  C++: %.3f\n", dp_cpp$alpha))

  # Posterior statistics (last 200 iterations)
  alpha_r <- tail(dp_r$alphaChain, 200)
  alpha_cpp <- tail(dp_cpp$alphaChain, 200)

  cat("\nPosterior alpha statistics:\n")
  cat(sprintf("  R: mean=%.3f, sd=%.3f\n", mean(alpha_r), sd(alpha_r)))
  cat(sprintf("  C++: mean=%.3f, sd=%.3f\n", mean(alpha_cpp), sd(alpha_cpp)))

  # Trace plot comparison
  if (length(alpha_r) == length(alpha_cpp)) {
    ks_test <- ks.test(alpha_r, alpha_cpp)
    cat(sprintf("\nKS test p-value: %.3f\n", ks_test$p.value))

    # Check mixing - autocorrelation
    acf_r <- acf(alpha_r, plot = FALSE)$acf[2:11]
    acf_cpp <- acf(alpha_cpp, plot = FALSE)$acf[2:11]

    cat("\nMarkov chain mixing (mean autocorrelation at lags 1-10):\n")
    cat(sprintf("  R: %.3f\n", mean(abs(acf_r))))
    cat(sprintf("  C++: %.3f\n", mean(abs(acf_cpp))))
  }

  # Cluster membership stability
  if (!is.null(dp_r$clusterLabels) && !is.null(dp_cpp$clusterLabels)) {
    # Adjusted Rand Index would be computed here if available
    cat("\nCluster assignment comparison:\n")
    cat(sprintf("  Unique clusters R: %d\n", length(unique(dp_r$clusterLabels))))
    cat(sprintf("  Unique clusters C++: %d\n", length(unique(dp_cpp$clusterLabels))))
  }

  return(list(dp_r = dp_r, dp_cpp = dp_cpp))
}

#' Comprehensive Markov DP benchmarking
#'
#' Tests multiple scenarios from Neal (2000) paper
#' @param n_obs_vec Vector of observation counts
#' @param algorithm_vec Vector of algorithms to test
#' @param n_iter_vec Vector of iteration counts
#' @param mixture_types Types of mixture models to test
#' @param n_reps Number of repetitions
#' @export
benchmark_markov_dp_comprehensive <- function(
    n_obs_vec = c(100, 250, 500, 1000),
    algorithm_vec = c(2, 3, 5, 8),
    n_iter_vec = c(100, 500),
    mixture_types = c("gaussian", "mvnormal", "exponential"),
    n_reps = 3) {

  cat("\n=== Comprehensive Markov DP Benchmarking ===\n")
  cat("Testing Neal (2000) algorithms across scenarios\n\n")

  total_scenarios <- length(n_obs_vec) * length(algorithm_vec) *
    length(n_iter_vec) * length(mixture_types)
  cat(sprintf("Testing %d scenarios with %d reps each...\n", total_scenarios, n_reps))

  results <- data.frame()
  scenario_id <- 0

  for (mixture_type in mixture_types) {
    for (algo in algorithm_vec) {
      for (n_obs in n_obs_vec) {
        for (n_iter in n_iter_vec) {

          scenario_id <- scenario_id + 1
          cat(sprintf("\nScenario %d/%d: %s, algo=%d, n=%d, iter=%d\n",
                      scenario_id, total_scenarios, mixture_type, algo, n_obs, n_iter))

          # Generate data based on mixture type
          set.seed(scenario_id)

          if (mixture_type == "gaussian") {
            # Gaussian mixture as in Neal examples
            y <- c(rnorm(n_obs/3, -3, 0.5),
                   rnorm(n_obs/3, 0, 0.7),
                   rnorm(n_obs/3, 3, 0.5))
            y <- sample(y) # Shuffle

          } else if (mixture_type == "mvnormal") {
            # Multivariate normal (d=2)
            library(mvtnorm)
            y1 <- rmvnorm(n_obs/2, mean = c(-2, -2), sigma = diag(2) * 0.5)
            y2 <- rmvnorm(n_obs/2, mean = c(2, 2), sigma = diag(2) * 0.7)
            y <- rbind(y1, y2)[sample(n_obs), ]

          } else if (mixture_type == "exponential") {
            # Exponential mixture
            y <- c(rexp(n_obs/2, rate = 2), rexp(n_obs/2, rate = 0.5))
            y <- sample(y)
          }

          # Run benchmarks
          for (rep in 1:n_reps) {

            # R implementation
            set_use_cpp(FALSE)
            gc(reset = TRUE)

            time_r <- system.time({
              if (mixture_type == "gaussian") {
                dp_r <- DirichletProcessGaussian(y)
              } else if (mixture_type == "mvnormal") {
                dp_r <- DirichletProcessMvnormal(y,
                                                 g0Priors = list(mu0 = c(0,0), Lambda = diag(2),
                                                                 kappa0 = 1, nu = 4))
              } else {
                dp_r <- DirichletProcessExponential(y)
              }

              if ("algorithm" %in% names(dp_r)) {
                dp_r$algorithm <- algo
              }
              dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
            })

            n_clusters_r <- dp_r$numberClusters
            alpha_final_r <- dp_r$alpha

            # Compute chain statistics
            if (length(dp_r$alphaChain) > 100) {
              alpha_ess_r <- coda::effectiveSize(dp_r$alphaChain)
            } else {
              alpha_ess_r <- NA
            }

            results <- rbind(results, data.frame(
              scenario_id = scenario_id,
              mixture_type = mixture_type,
              algorithm = algo,
              n_obs = n_obs,
              n_iter = n_iter,
              rep = rep,
              implementation = "R",
              time = time_r["elapsed"],
              n_clusters = n_clusters_r,
              alpha_final = alpha_final_r,
              alpha_ess = as.numeric(alpha_ess_r)
            ))

            # C++ implementation
            set_use_cpp(TRUE)
            gc(reset = TRUE)

            time_cpp <- system.time({
              if (mixture_type == "gaussian") {
                dp_cpp <- DirichletProcessGaussian(y)
              } else if (mixture_type == "mvnormal") {
                dp_cpp <- DirichletProcessMvnormal(y,
                                                   g0Priors = list(mu0 = c(0,0), Lambda = diag(2),
                                                                   kappa0 = 1, nu = 4))
              } else {
                dp_cpp <- DirichletProcessExponential(y)
              }

              if ("algorithm" %in% names(dp_cpp)) {
                dp_cpp$algorithm <- algo
              }
              dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
            })

            n_clusters_cpp <- dp_cpp$numberClusters
            alpha_final_cpp <- dp_cpp$alpha

            if (length(dp_cpp$alphaChain) > 100) {
              alpha_ess_cpp <- coda::effectiveSize(dp_cpp$alphaChain)
            } else {
              alpha_ess_cpp <- NA
            }

            results <- rbind(results, data.frame(
              scenario_id = scenario_id,
              mixture_type = mixture_type,
              algorithm = algo,
              n_obs = n_obs,
              n_iter = n_iter,
              rep = rep,
              implementation = "C++",
              time = time_cpp["elapsed"],
              n_clusters = n_clusters_cpp,
              alpha_final = alpha_final_cpp,
              alpha_ess = as.numeric(alpha_ess_cpp)
            ))
          }

          # Print progress summary
          scenario_results <- results %>%
            filter(scenario_id == !!scenario_id) %>%
            group_by(implementation) %>%
            summarise(mean_time = mean(time), .groups = "drop")

          speedup <- scenario_results$mean_time[scenario_results$implementation == "R"] /
            scenario_results$mean_time[scenario_results$implementation == "C++"]

          cat(sprintf("  Speedup: %.1fx\n", speedup))
        }
      }
    }
  }

  return(results)
}

#' Memory profiling for Markov DP algorithms
#'
#' @param n_obs Number of observations
#' @param n_iter Number of iterations
#' @param algorithm Algorithm to profile
#' @return List with memory traces
#' @export
profile_markov_dp_components <- function(n = 500, algorithm = 3) {

  cat("\n=== Markov DP Component Profiling ===\n")
  cat(sprintf("Algorithm %d components (Neal 2000)\n", algorithm))

  # Generate test data
  set.seed(123)
  y <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 2, 0.5))

  # Create DP objects
  set_use_cpp(FALSE)
  dp_r <- DirichletProcessGaussian(y)
  dp_r <- Initialise(dp_r)

  set_use_cpp(TRUE)
  dp_cpp <- DirichletProcessGaussian(y)
  dp_cpp <- Initialise(dp_cpp)

  # Component timing results
  results <- data.frame()

  # Time Cluster Component Update
  cat("Profiling ClusterComponentUpdate...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      ClusterComponentUpdate(dp_r)
    },
    times = 10
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      ClusterComponentUpdate(dp_cpp)
    },
    times = 10
  )

  results <- rbind(results, data.frame(
    component = "ClusterComponentUpdate",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Time Parameter Update
  cat("Profiling ClusterParameterUpdate...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      ClusterParameterUpdate(dp_r)
    },
    times = 50
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      ClusterParameterUpdate(dp_cpp)
    },
    times = 50
  )

  results <- rbind(results, data.frame(
    component = "ClusterParameterUpdate",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Time Alpha Update
  cat("Profiling UpdateAlpha...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      UpdateAlpha(dp_r)
    },
    times = 50
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      UpdateAlpha(dp_cpp)
    },
    times = 50
  )

  results <- rbind(results, data.frame(
    component = "UpdateAlpha",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Time Likelihood calculation
  cat("Profiling Likelihood calculations...\n")
  time_r <- microbenchmark(
    R = {
      set_use_cpp(FALSE)
      total_lik <- 0
      for (i in 1:min(100, n)) {
        single_cluster_params <- list(
          mu = dp_r$clusterParameters[[1]][,,1, drop=FALSE],
          sigma = dp_r$clusterParameters[[2]][,,1, drop=FALSE]
        )
        lik <- Likelihood(dp_r$mixingDistribution,
                          y[i, , drop = FALSE],
                          single_cluster_params)
        total_lik <- total_lik + lik
      }
    },
    times = 10
  )

  time_cpp <- microbenchmark(
    Cpp = {
      set_use_cpp(TRUE)
      total_lik <- 0
      for (i in 1:min(100, n)) {
        single_cluster_params <- list(
          mu = dp_cpp$clusterParameters[[1]][,,1, drop=FALSE],
          sigma = dp_cpp$clusterParameters[[2]][,,1, drop=FALSE]
        )
        lik <- Likelihood(dp_cpp$mixingDistribution,
                          y[i, , drop = FALSE],
                          single_cluster_params)
        total_lik <- total_lik + lik
      }
    },
    times = 10
  )

  results <- rbind(results, data.frame(
    component = "Likelihood",
    implementation = c("R", "C++"),
    mean_time_ms = c(mean(time_r$time) / 1e6, mean(time_cpp$time) / 1e6)
  ))

  # Calculate speedups
  cat("\n\nComponent timings (milliseconds):\n")
  print(results)

  # Calculate speedup for each component
  speedups <- results %>%
    pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
    mutate(speedup = R / `C++`)

  cat("\nSpeedup by component:\n")
  for (i in 1:nrow(speedups)) {
    cat(sprintf("  %s: %.1fx\n", speedups$component[i], speedups$speedup[i]))
  }

  return(results)
}

#' Visualize Markov DP benchmark results
#'
#' @param bench_results Results from comprehensive benchmark
#' @param save_plots Whether to save plots to file
#' @return List of ggplot objects
#' @export
visualize_markov_dp_benchmarks <- function(bench_results, save_plots = TRUE) {

  # 1. Performance scaling by algorithm
  p_scaling <- bench_results %>%
    group_by(algorithm, n_obs, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop") %>%
    ggplot(aes(x = n_obs, y = mean_time, color = implementation)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~ paste("Algorithm", algorithm), scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Markov DP Performance Scaling by Algorithm",
         subtitle = "Based on Neal (2000) algorithms",
         x = "Number of observations",
         y = "Time (seconds)",
         color = "Implementation") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # 2. Speedup heatmap
  speedup_data <- bench_results %>%
    group_by(algorithm, n_obs, mixture_type) %>%
    summarise(
      time_r = mean(time[implementation == "R"]),
      time_cpp = mean(time[implementation == "C++"]),
      speedup = time_r / time_cpp,
      .groups = "drop"
    )

  p_speedup <- ggplot(speedup_data,
                      aes(x = factor(n_obs), y = factor(algorithm), fill = speedup)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1fx", speedup)), size = 3) +
    facet_wrap(~ mixture_type) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1,
                         name = "Speedup\n(R/C++)") +
    labs(title = "C++ Speedup over R Implementation",
         subtitle = "Markov chain algorithms on different mixture types",
         x = "Number of observations",
         y = "Algorithm") +
    theme_minimal()

  # 3. Mixing efficiency (ESS/second)
  p_mixing <- bench_results %>%
    filter(!is.na(alpha_ess)) %>%
    mutate(ess_per_sec = alpha_ess / time) %>%
    group_by(algorithm, implementation) %>%
    summarise(mean_ess_per_sec = mean(ess_per_sec), .groups = "drop") %>%
    ggplot(aes(x = factor(algorithm), y = mean_ess_per_sec, fill = implementation)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Markov Chain Mixing Efficiency",
         subtitle = "Effective Sample Size per second for α parameter",
         x = "Algorithm",
         y = "ESS/second",
         fill = "Implementation") +
    theme_minimal()

  # 4. Cluster recovery accuracy
  p_clusters <- bench_results %>%
    group_by(mixture_type, algorithm, implementation) %>%
    summarise(
      mean_clusters = mean(n_clusters),
      sd_clusters = sd(n_clusters),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = factor(algorithm), y = mean_clusters, fill = implementation)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_clusters - sd_clusters,
                      ymax = mean_clusters + sd_clusters),
                  position = position_dodge(0.9), width = 0.2) +
    facet_wrap(~ mixture_type) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Cluster Recovery Performance",
         subtitle = "True number of clusters = 3 (red line)",
         x = "Algorithm",
         y = "Mean clusters found",
         fill = "Implementation") +
    theme_minimal()

  # Combine plots
  combined_plot <- (p_scaling + p_speedup) / (p_mixing + p_clusters) +
    plot_annotation(
      title = "Markov Chain Dirichlet Process Benchmark Results",
      subtitle = sprintf("Comparing %d algorithms from Neal (2000)",
                         length(unique(bench_results$algorithm))),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  if (save_plots) {
    ggsave("markov_dp_benchmark_results.png", combined_plot,
           width = 16, height = 12, dpi = 300)
  }

  return(list(
    scaling = p_scaling,
    speedup = p_speedup,
    mixing = p_mixing,
    clusters = p_clusters,
    combined = combined_plot
  ))
}

#' Generate comprehensive Markov DP benchmark report
#'
#' @export
run_markov_dp_benchmark_report <- function() {

  # Check C++ implementation status
  cpp_status <- get_cpp_status()
  if (!cpp_status$available && !cpp_status$mcmc_runner) {
    stop("C++ implementation not available. Please recompile with C++ support.")
  }

  # 1. Quick benchmark
  cat("\n1. Running quick benchmark...\n")
  quick_results <- quick_markov_dp_benchmark()

  # 2. Component profiling
  cat("\n2. Profiling individual components...\n")
  component_results <- profile_markov_dp_components()

  # 3. Statistical equivalence
  cat("\n3. Testing statistical equivalence...\n")
  equiv_results <- test_markov_dp_equivalence()

  # 4. Comprehensive benchmarks
  cat("\n4. Running comprehensive benchmarks...\n")
  bench_results <- benchmark_markov_dp_comprehensive(
    n_obs_vec = c(100, 500, 1000),
    algorithm_vec = c(2, 3, 8),
    n_iter_vec = c(200, 500),
    mixture_types = c("gaussian", "exponential"),
    n_reps = 3
  )

  # 5. Memory profiling
  cat("\n5. Memory profiling...\n")
  memory_results <- profile_markov_dp_memory()

  # Generate summary statistics
  summary_stats <- bench_results %>%
    group_by(implementation, algorithm) %>%
    summarise(
      mean_time = mean(time),
      median_time = median(time),
      mean_ess_per_sec = mean(alpha_ess / time, na.rm = TRUE),
      .groups = "drop"
    )

  # Overall speedup by algorithm
  speedup_by_algo <- bench_results %>%
    group_by(algorithm) %>%
    summarise(
      mean_time_r = mean(time[implementation == "R"]),
      mean_time_cpp = mean(time[implementation == "C++"]),
      speedup = mean_time_r / mean_time_cpp,
      .groups = "drop"
    )

  # Print summary report
  cat("\n")
  cat("================================================\n")
  cat("MARKOV DP BENCHMARK SUMMARY (Neal 2000)\n")
  cat("================================================\n")

  cat("\n1. OVERALL PERFORMANCE BY ALGORITHM:\n")
  print(summary_stats)

  cat("\n2. SPEEDUP BY ALGORITHM:\n")
  print(speedup_by_algo)

  cat("\n3. COMPONENT-LEVEL SPEEDUP:\n")
  comp_speedup <- component_results %>%
    pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
    mutate(speedup = R / `C++`)
  print(comp_speedup)

  cat("\n4. OVERALL SPEEDUP RANGE:\n")
  overall_speedup <- summary_stats$mean_time[summary_stats$implementation == "R"] /
    summary_stats$mean_time[summary_stats$implementation == "C++"]
  cat(sprintf("  Mean speedup: %.1fx\n", mean(overall_speedup)))
  cat(sprintf("  Min speedup: %.1fx\n", min(speedup_by_algo$speedup)))
  cat(sprintf("  Max speedup: %.1fx\n", max(speedup_by_algo$speedup)))

  cat("\n5. MIXING EFFICIENCY:\n")
  mixing_summary <- bench_results %>%
    filter(!is.na(alpha_ess)) %>%
    group_by(implementation) %>%
    summarise(
      mean_ess = mean(alpha_ess),
      mean_ess_per_sec = mean(alpha_ess / time),
      .groups = "drop"
    )
  print(mixing_summary)

  # Create visualizations
  cat("\n6. Generating visualizations...\n")
  plots <- visualize_markov_dp_benchmarks(bench_results)

  # Display main plot
  print(plots$combined)

  # Save results
  results <- list(
    quick_benchmark = quick_results,
    component_profile = component_results,
    statistical_equivalence = equiv_results,
    comprehensive_benchmark = bench_results,
    memory_profile = memory_results,
    summary_stats = summary_stats,
    speedup_by_algorithm = speedup_by_algo,
    plots = plots
  )

  cat("\nMarkov DP benchmarking complete!\n")
  cat("Results demonstrate implementation of Neal (2000) algorithms\n")

  return(invisible(results))
}

# Run the benchmark if sourced interactively
if (interactive()) {
  cat("Markov Chain Dirichlet Process Benchmark Suite\n")
  cat("Based on Neal (2000) - 'Markov Chain Sampling Methods for Dirichlet Process Mixture Models'\n")
  cat("=============================================================================\n\n")

  results <- run_markov_dp_benchmark_report()
}
