#' Benchmark Hierarchical MVNormal implementations
#'
#' This script compares the performance of R and C++ implementations
#' of hierarchical Dirichlet Process with MVNormal base distributions

library(dirichletprocess)
library(ggplot2)
library(microbenchmark)
library(mvtnorm)
library(dplyr)
library(gridExtra)

#' Quick benchmark for hierarchical MVNormal
#'
#' @param n_groups Vector of group counts to test
#' @param n_obs_per_group Vector of observations per group
#' @param dimensions Vector of data dimensions
#' @param n_iter Number of MCMC iterations
#' @param n_reps Number of repetitions
#' @export
quick_hierarchical_mvnormal_benchmark <- function(
    n_groups = c(2, 3),
    n_obs_per_group = c(30, 50),
    dimensions = c(2, 3),
    n_iter = 30,
    n_reps = 2) {

  cat("=== Quick Hierarchical MVNormal Distribution Benchmark ===\n")

  results <- data.frame()

  for (g in n_groups) {
    for (n in n_obs_per_group) {
      for (d in dimensions) {
        cat(sprintf("Testing groups = %d, n/group = %d, d = %d:\n", g, n, d))

        # Generate test data with hierarchical structure
        set.seed(123 + g * 100 + n * 10 + d)
        dataList <- list()

        for (i in 1:g) {
          # Each group has its own mean
          mu_group <- rep((i-1) * 3, d)
          dataList[[i]] <- rmvnorm(n, mu_group, diag(d))
        }

        # Benchmark R implementation
        times_r <- numeric(n_reps)
        gamma_final_r <- numeric(n_reps)

        for (i in 1:n_reps) {
          # Set up priors for R implementation
          g0Priors_R <- list(
            mu0 = rep(0, d),
            phi0 = diag(d),
            sigma0 = diag(d),
            nu0 = d + 2
          )

          if (exists("enable_cpp_hierarchical_samplers")) {
            enable_cpp_hierarchical_samplers(FALSE)
          }
          set_use_cpp(FALSE)
          cat("Using R implementations for hierarchical samplers\n")

          time_r <- system.time({
            hdp_r <- DirichletProcessHierarchicalMvnormal2(
              dataList = dataList,
              g0Priors = g0Priors_R,
              gammaPriors = c(2, 4),
              alphaPriors = c(2, 4),
              numSticks = 20,
              numInitialClusters = 1
            )

            # Fix the class order to ensure hierarchical dispatch
            class(hdp_r) <- c("hierarchical", "dirichletprocess", "list")

            hdp_r <- Fit(hdp_r, n_iter, progressBar = FALSE)
          })
          times_r[i] <- time_r["elapsed"]
          gamma_final_r[i] <- hdp_r$gamma
        }

        # Benchmark C++ implementation (if available)
        times_cpp <- numeric(n_reps)
        gamma_final <- numeric(n_reps)

        cpp_available <- exists("hierarchical_mvnormal_run")

        if (cpp_available) {
          for (i in 1:n_reps) {
            # Set up priors for C++ implementation
            g0Priors_CPP <- list(
              mu0 = rep(0, d),
              kappa0 = 1.0,
              Lambda = diag(d),
              nu = d + 2
            )

            if (exists("enable_cpp_hierarchical_samplers")) {
              enable_cpp_hierarchical_samplers(TRUE)
            }
            set_use_cpp(TRUE)
            cat("C++ samplers enabled for hierarchical Dirichlet processes\n")

            time_cpp <- system.time({
              hdp_params <- list(
                n_sticks = 20,
                prior_params = g0Priors_CPP,
                alpha_prior = c(2, 4),
                gamma_prior = c(2, 4)
              )

              mcmc_params <- list(
                n_iter = n_iter,
                n_burn = 0,
                thin = 1,
                update_prior = TRUE,
                show_progress = FALSE
              )

              result <- hierarchical_mvnormal_run(dataList, hdp_params, mcmc_params)
            })
            times_cpp[i] <- time_cpp["elapsed"]
            gamma_final[i] <- result$final_state$gamma
          }
        } else {
          cat("  C++ implementation not available\n")
          times_cpp <- rep(NA, n_reps)
          gamma_final <- rep(NA, n_reps)
        }

        # Calculate speedup
        speedup <- ifelse(!is.na(times_cpp[1]),
                          mean(times_r) / mean(times_cpp),
                          NA)

        cat(sprintf("  R:   %.3f sec (gamma = %.3f)\n",
                    mean(times_r), mean(gamma_final_r)))
        if (!is.na(speedup)) {
          cat(sprintf("  C++: %.3f sec (gamma = %.3f)\n",
                      mean(times_cpp), mean(gamma_final)))
          cat(sprintf("  Speedup: %.1fx\n", speedup))
        }

        # Store results
        results <- rbind(results, data.frame(
          n_groups = g,
          n_per_group = n,
          dimension = d,
          implementation = "R",
          time = times_r,
          memory_mb = NA,
          gamma = gamma_final_r
        ))

        results <- rbind(results, data.frame(
          n_groups = g,
          n_per_group = n,
          dimension = d,
          implementation = "C++",
          time = times_cpp,
          memory_mb = NA,
          gamma = gamma_final
        ))
      }
    }
  }

  return(results)
}

#' Comprehensive Hierarchical MVNormal benchmarking
#'
#' @param n_groups_vec Vector of group counts
#' @param n_obs_vec Vector of observations per group
#' @param dim_vec Vector of dimensions
#' @param n_iter_vec Vector of iteration counts
#' @param n_reps Number of repetitions
#' @export
benchmark_hierarchical_mvnormal_comprehensive <- function(
    n_groups_vec = c(2, 3, 5),
    n_obs_vec = c(30, 50, 100),
    dim_vec = c(2, 3),
    n_iter_vec = c(50, 100),
    n_reps = 3) {

  cat("\n=== Comprehensive Hierarchical MVNormal Benchmarking ===\n\n")

  total_scenarios <- length(n_groups_vec) * length(n_obs_vec) *
    length(dim_vec) * length(n_iter_vec)
  cat(sprintf("Testing %d scenarios with %d reps each...\n", total_scenarios, n_reps))

  results <- data.frame()
  scenario_id <- 0

  for (n_groups in n_groups_vec) {
    for (n_obs in n_obs_vec) {
      for (d in dim_vec) {
        for (n_iter in n_iter_vec) {

          scenario_id <- scenario_id + 1
          cat(sprintf("\nScenario %d/%d: groups=%d, n/group=%d, d=%d, iter=%d\n",
                      scenario_id, total_scenarios, n_groups, n_obs, d, n_iter))

          # Generate hierarchical data
          set.seed(scenario_id)
          dataList <- list()

          for (g in 1:n_groups) {
            # Create 2-3 clusters per group
            n_clusters_g <- sample(2:3, 1)
            data_g <- matrix(0, n_obs, d)
            n_per_cluster <- n_obs / n_clusters_g

            for (k in 1:n_clusters_g) {
              start_idx <- floor((k-1) * n_per_cluster) + 1
              end_idx <- min(floor(k * n_per_cluster), n_obs)

              mu_k <- rep((g-1) * 5 + k * 3, d)
              data_g[start_idx:end_idx, ] <- rmvnorm(
                end_idx - start_idx + 1,
                mu_k,
                diag(d) * 0.5
              )
            }
            dataList[[g]] <- data_g
          }

          # Benchmark both implementations
          for (impl in c("R", "C++")) {

            times <- numeric(n_reps)
            gamma_values <- numeric(n_reps)

            for (rep in 1:n_reps) {

              if (impl == "R") {
                # Set up priors for R implementation
                g0Priors_R <- list(
                  mu0 = rep(0, d),
                  phi0 = diag(d),
                  sigma0 = diag(d),
                  nu0 = d + 2
                )

                set_use_cpp(FALSE)
                time_taken <- system.time({
                  hdp <- DirichletProcessHierarchicalMvnormal2(
                    dataList = dataList,
                    g0Priors = g0Priors_R,
                    gammaPriors = c(2, 4),
                    alphaPriors = c(2, 4),
                    numSticks = 30
                  )
                  hdp <- Fit(hdp, n_iter, progressBar = FALSE)
                })
                times[rep] <- time_taken["elapsed"]
                gamma_values[rep] <- hdp$gamma

              } else {
                # Set up priors for C++ implementation
                g0Priors_CPP <- list(
                  mu0 = rep(0, d),
                  kappa0 = 0.5,
                  Lambda = diag(d),
                  nu = d + 2
                )

                set_use_cpp(TRUE)
                hdp_params <- list(
                  n_sticks = 30,
                  prior_params = g0Priors_CPP,
                  alpha_prior = c(2, 4),
                  gamma_prior = c(2, 4)
                )

                mcmc_params <- list(
                  n_iter = n_iter,
                  n_burn = 0,
                  thin = 1,
                  update_prior = TRUE,
                  show_progress = FALSE
                )

                time_taken <- system.time({
                  result <- hierarchical_mvnormal_run(dataList, hdp_params, mcmc_params)
                })
                times[rep] <- time_taken["elapsed"]
                gamma_values[rep] <- result$final_state$gamma
              }
            }

            # Store results
            results <- rbind(results, data.frame(
              n_groups = n_groups,
              n_obs_per_group = n_obs,
              dimension = d,
              n_iter = n_iter,
              implementation = impl,
              time = times,
              gamma = gamma_values,
              memory_mb = NA
            ))
          }
        }
      }
    }
  }

  return(results)
}

#' Profile Hierarchical MVNormal components
#'
#' @param n_groups Number of groups
#' @param n_per_group Observations per group
#' @param d Data dimension
#' @export
profile_hierarchical_mvnormal_components <- function(n_groups = 3, n_per_group = 50, d = 3) {

  cat("\n=== Hierarchical MVNormal Component Profiling ===\n\n")

  # Generate test data
  set.seed(123)
  dataList <- list()
  for (i in 1:n_groups) {
    mu_g <- rep((i-1) * 3, d)
    dataList[[i]] <- rmvnorm(n_per_group, mu_g, diag(d))
  }

  # Set up priors for R implementation
  g0Priors_R <- list(
    mu0 = rep(0, d),
    phi0 = diag(d),
    sigma0 = diag(d),
    nu0 = d + 2
  )

  # Create hierarchical DP
  hdp <- DirichletProcessHierarchicalMvnormal2(
    dataList = dataList,
    g0Priors = g0Priors_R,
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    numSticks = 20
  )

  # Fix class order
  class(hdp) <- c("hierarchical", "dirichletprocess", "list")

  # Profile individual components
  cat("Profiling component updates:\n")

  # 1. Cluster component update
  time_comp <- system.time({
    for (i in 1:10) {
      hdp <- ClusterComponentUpdate(hdp)
    }
  })
  cat(sprintf("  Cluster component update: %.3f ms/iter\n", time_comp["elapsed"] * 100))

  # 2. Alpha update
  time_alpha <- system.time({
    for (i in 1:10) {
      hdp <- UpdateAlpha(hdp)
    }
  })
  cat(sprintf("  Alpha update: %.3f ms/iter\n", time_alpha["elapsed"] * 100))

  # 3. Global parameter update
  time_global <- system.time({
    for (i in 1:10) {
      hdp <- GlobalParameterUpdate(hdp)
    }
  })
  cat(sprintf("  Global parameter update: %.3f ms/iter\n", time_global["elapsed"] * 100))

  # 4. Gamma update
  time_gamma <- system.time({
    for (i in 1:10) {
      hdp <- UpdateGamma(hdp)
    }
  })
  cat(sprintf("  Gamma update: %.3f ms/iter\n", time_gamma["elapsed"] * 100))

  # Profile full iteration
  time_full <- system.time({
    hdp_copy <- hdp
    class(hdp_copy) <- c("hierarchical", "dirichletprocess", "list")
    hdp_copy <- Fit(hdp_copy, 10, progressBar = FALSE)
  })
  cat(sprintf("\nFull iteration: %.3f ms/iter\n", time_full["elapsed"] * 100))

  return(list(
    component_update = time_comp["elapsed"] / 10,
    alpha_update = time_alpha["elapsed"] / 10,
    global_update = time_global["elapsed"] / 10,
    gamma_update = time_gamma["elapsed"] / 10,
    full_iteration = time_full["elapsed"] / 10
  ))
}

#' Test statistical equivalence between R and C++ implementations
#'
#' @param n_groups Number of groups
#' @param n_per_group Observations per group
#' @param d Dimension
#' @param n_iter Number of iterations
#' @export
test_hierarchical_mvnormal_statistical_equivalence <- function(n_groups = 2,
                                                               n_per_group = 30,
                                                               d = 2,
                                                               n_iter = 100) {

  cat("\n=== Testing Statistical Equivalence ===\n\n")

  # Generate test data
  set.seed(456)
  dataList <- list()

  for (g in 1:n_groups) {
    mu_g <- rep(g * 2, d)
    dataList[[g]] <- rmvnorm(n_per_group, mu_g, diag(d))
  }

  # Run R implementation
  set.seed(789)
  set_use_cpp(FALSE)

  # Set up priors for R implementation
  g0Priors_R <- list(
    mu0 = rep(0, d),
    phi0 = diag(d),
    sigma0 = diag(d),
    nu0 = d + 2
  )

  hdp_r <- DirichletProcessHierarchicalMvnormal2(
    dataList = dataList,
    g0Priors = g0Priors_R
  )
  hdp_r <- Fit(hdp_r, n_iter, progressBar = FALSE)

  # Run C++ implementation
  set.seed(789)
  set_use_cpp(TRUE)

  # Set up priors for C++ implementation
  g0Priors_CPP <- list(
    mu0 = rep(0, d),
    kappa0 = 1.0,
    Lambda = diag(d),
    nu = d + 2
  )

  hdp_params <- list(
    n_sticks = 50,
    prior_params = g0Priors_CPP,
    alpha_prior = c(2, 4),
    gamma_prior = c(2, 4)
  )

  mcmc_params <- list(
    n_iter = n_iter,
    n_burn = 0,
    thin = 1,
    update_prior = TRUE,
    show_progress = FALSE
  )

  result_cpp <- hierarchical_mvnormal_run(dataList, hdp_params, mcmc_params)

  # Compare results
  cat("Number of clusters per group:\n")
  cat("  R implementation:\n")
  for (i in 1:n_groups) {
    cat(sprintf("    Group %d: %d clusters\n", i, hdp_r$indDP[[i]]$numberClusters))
  }

  cat("\n  C++ implementation:\n")
  final_labels <- result_cpp$final_state$cluster_labels
  for (i in 1:n_groups) {
    n_clusters_cpp <- length(unique(final_labels[[i]]))
    cat(sprintf("    Group %d: %d clusters\n", i, n_clusters_cpp))
  }

  # Compare gamma values
  cat("\nGlobal concentration parameter (gamma):\n")
  cat(sprintf("  R: %.3f\n", hdp_r$gamma))
  cat(sprintf("  C++: %.3f\n", result_cpp$final_state$gamma))

  # Compare alpha values
  cat("\nLocal concentration parameters (alpha):\n")
  cat("  R:", sapply(hdp_r$indDP, function(x) round(x$alpha, 3)), "\n")
  cat("  C++:", round(result_cpp$final_state$alphas, 3), "\n")

  return(list(hdp_r = hdp_r, result_cpp = result_cpp))
}

#' Visualize Hierarchical MVNormal benchmark results
#'
#' @param bench_results Results from comprehensive benchmarking
#' @export
visualize_hierarchical_mvnormal_benchmarks <- function(bench_results) {

  # Implementation colors
  impl_colors <- c("R" = "#E41A1C", "C++" = "#377EB8")

  # 1. Speedup heatmap
  speedup_data <- bench_results %>%
    group_by(n_groups, n_obs_per_group, dimension, n_iter) %>%
    summarise(
      time_r = mean(time[implementation == "R"]),
      time_cpp = mean(time[implementation == "C++"]),
      speedup = time_r / time_cpp,
      .groups = "drop"
    )

  p_speedup <- ggplot(speedup_data,
                      aes(x = factor(n_obs_per_group),
                          y = factor(n_groups),
                          fill = speedup)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1fx", speedup)), size = 3) +
    facet_grid(dimension ~ n_iter,
               labeller = labeller(dimension = label_both, n_iter = label_both)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 1, limits = c(0, NA)) +
    labs(
      title = "C++ Speedup over R Implementation",
      x = "Observations per Group",
      y = "Number of Groups",
      fill = "Speedup"
    ) +
    theme_minimal()

  # 2. Time comparison by scenario
  p_time <- ggplot(bench_results,
                   aes(x = factor(n_obs_per_group), y = time,
                       fill = implementation)) +
    geom_boxplot(alpha = 0.8) +
    facet_grid(n_groups ~ dimension,
               labeller = labeller(n_groups = label_both, dimension = label_both)) +
    scale_fill_manual(values = impl_colors) +
    labs(
      title = "Execution Time by Scenario",
      x = "Observations per Group",
      y = "Time (seconds)",
      fill = "Implementation"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # 3. Scaling analysis
  scaling_data <- bench_results %>%
    filter(dimension == 2, n_iter == min(n_iter)) %>%
    group_by(n_groups, n_obs_per_group, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop")

  p_scaling <- ggplot(scaling_data,
                      aes(x = n_obs_per_group, y = mean_time,
                          color = implementation, linetype = factor(n_groups))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Scaling with Data Size",
      x = "Observations per Group",
      y = "Mean Time (seconds)",
      color = "Implementation",
      linetype = "Groups"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # 4. Convergence comparison (gamma parameter)
  gamma_data <- bench_results %>%
    group_by(n_groups, dimension, implementation) %>%
    summarise(
      mean_gamma = mean(gamma),
      sd_gamma = sd(gamma),
      .groups = "drop"
    )

  p_gamma <- ggplot(gamma_data,
                    aes(x = factor(n_groups), y = mean_gamma,
                        fill = implementation)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean_gamma - sd_gamma,
                      ymax = mean_gamma + sd_gamma),
                  position = position_dodge(0.9), width = 0.2) +
    facet_wrap(~dimension, labeller = labeller(dimension = label_both)) +
    scale_fill_manual(values = impl_colors) +
    labs(
      title = "Global Concentration Parameter (γ) Estimates",
      x = "Number of Groups",
      y = "Mean γ",
      fill = "Implementation"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Create dashboard
  dashboard <- grid.arrange(
    p_speedup, p_time,
    p_scaling, p_gamma,
    ncol = 2,
    top = "Hierarchical MVNormal Performance Comparison"
  )

  return(list(
    speedup = p_speedup,
    time = p_time,
    scaling = p_scaling,
    gamma = p_gamma,
    dashboard = dashboard
  ))
}

#' Run complete hierarchical MVNormal benchmark report
#'
#' @export
run_hierarchical_mvnormal_benchmark_report <- function() {

  cat("\n")
  cat("================================================\n")
  cat("Hierarchical MVNormal Dirichlet Process Benchmarking\n")
  cat("================================================\n")

  # 1. Quick benchmark
  cat("\n1. Running quick benchmark...\n")
  quick_results <- quick_hierarchical_mvnormal_benchmark(
    n_groups = c(2, 3),
    n_obs_per_group = c(30, 50),
    dimensions = c(2, 3),
    n_iter = 30,
    n_reps = 2
  )

  # 2. Component profiling
  cat("\n2. Profiling individual components...\n")
  component_results <- profile_hierarchical_mvnormal_components()

  # 3. Statistical equivalence
  cat("\n3. Testing statistical equivalence...\n")
  equiv_results <- test_hierarchical_mvnormal_statistical_equivalence(
    n_groups = 2,
    n_per_group = 30,
    d = 2,
    n_iter = 50
  )

  # 4. Comprehensive benchmarks
  cat("\n4. Running comprehensive benchmarks...\n")
  bench_results <- benchmark_hierarchical_mvnormal_comprehensive(
    n_groups_vec = c(2, 3),
    n_obs_vec = c(30, 50),
    dim_vec = c(2, 3),
    n_iter_vec = c(50),
    n_reps = 2
  )

  # Generate summary statistics
  summary_stats <- bench_results %>%
    group_by(implementation) %>%
    summarise(
      mean_time = mean(time),
      median_time = median(time),
      mean_gamma = mean(gamma),
      .groups = "drop"
    )

  # Calculate speedup by configuration
  speedup_by_config <- bench_results %>%
    group_by(n_groups, dimension) %>%
    summarise(
      mean_time_r = mean(time[implementation == "R"]),
      mean_time_cpp = mean(time[implementation == "C++"]),
      speedup = mean_time_r / mean_time_cpp,
      .groups = "drop"
    )

  # Print summary report
  cat("\n")
  cat("================================================\n")
  cat("HIERARCHICAL MVNORMAL BENCHMARK SUMMARY\n")
  cat("================================================\n")

  cat("\n1. OVERALL PERFORMANCE:\n")
  print(summary_stats)

  cat("\n2. SPEEDUP BY CONFIGURATION:\n")
  print(speedup_by_config)

  cat("\n3. OVERALL SPEEDUP:\n")
  overall_speedup <- summary_stats$mean_time[summary_stats$implementation == "R"] /
    summary_stats$mean_time[summary_stats$implementation == "C++"]
  cat(sprintf("  Mean speedup: %.1fx\n", overall_speedup))
  cat(sprintf("  Min speedup: %.1fx\n", min(speedup_by_config$speedup)))
  cat(sprintf("  Max speedup: %.1fx\n", max(speedup_by_config$speedup)))

  # Create visualizations
  cat("\n5. Generating visualizations...\n")
  plots <- visualize_hierarchical_mvnormal_benchmarks(bench_results)

  # Display main dashboard
  print(plots$dashboard)

  # Save results
  results <- list(
    quick_benchmark = quick_results,
    component_profile = component_results,
    statistical_equivalence = equiv_results,
    comprehensive_benchmark = bench_results,
    summary_stats = summary_stats,
    speedup_by_config = speedup_by_config,
    plots = plots
  )

  cat("\nBenchmarking complete!\n")

  return(invisible(results))
}

# Run the benchmark if sourced interactively
if (interactive()) {
  results <- run_hierarchical_mvnormal_benchmark_report()
}
