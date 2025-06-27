# benchmark/benchmark_hierarchical_beta.R
# Comprehensive benchmarking for Hierarchical Beta Distribution
# Compares R vs C++ implementations

library(dirichletprocess)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(microbenchmark)

# Add this helper function after the library statements
enable_cpp_hierarchical_samplers <- function(use_cpp = TRUE) {
  # Suppress warnings about fallback since we're explicitly controlling this
  suppressWarnings({
    if (use_cpp && exists("hierarchical_beta_fit_cpp")) {
      options(dirichletprocess.use_cpp_hierarchical = TRUE)
    } else {
      options(dirichletprocess.use_cpp_hierarchical = FALSE)
    }
  })
  invisible(use_cpp)
}

#' Generate synthetic hierarchical beta mixture data
#'
#' @param n_groups Number of groups/datasets
#' @param n_per_group Number of observations per group
#' @param shared_params Shared parameters across groups
#' @param unique_params List of unique parameters for each group
#' @param weights Mixing weights for each component
#' @return List of datasets
generate_hierarchical_beta_data <- function(n_groups = 3,
                                            n_per_group = 200,
                                            shared_params = list(c(0.25, 5)),
                                            unique_params = list(c(0.75, 6), c(0.4, 10)),
                                            weights = c(0.5, 0.5)) {

  dataList <- list()

  for (g in 1:n_groups) {
    # Each group has the shared component plus one unique component
    if (g <= length(unique_params)) {
      params <- list(shared_params[[1]], unique_params[[g]])
    } else {
      # If more groups than unique params, reuse some
      params <- list(shared_params[[1]], unique_params[[((g-1) %% length(unique_params)) + 1]])
    }

    # Generate data from mixture
    n_comp1 <- round(n_per_group * weights[1])
    n_comp2 <- n_per_group - n_comp1

    data_comp1 <- rbeta(n_comp1, params[[1]][1] * params[[1]][2],
                        (1 - params[[1]][1]) * params[[1]][2])
    data_comp2 <- rbeta(n_comp2, params[[2]][1] * params[[2]][2],
                        (1 - params[[2]][1]) * params[[2]][2])

    dataList[[g]] <- c(data_comp1, data_comp2)
  }

  return(dataList)
}

#' Quick Hierarchical Beta benchmark for interactive testing
#'
#' @param n_groups_vec Vector of number of groups to test
#' @param n_per_group_vec Vector of observations per group to test
#' @param n_iter Number of MCMC iterations
#' @param n_reps Number of repetitions for timing
#' @export
quick_hierarchical_beta_benchmark <- function(n_groups_vec = c(2, 3, 5),
                                              n_per_group_vec = c(100, 200),
                                              n_iter = 100,
                                              n_reps = 3) {

  cat("\n=== Quick Hierarchical Beta Distribution Benchmark ===\n\n")

  results <- data.frame()

  for (n_groups in n_groups_vec) {
    for (n_per_group in n_per_group_vec) {
      cat(sprintf("Testing %d groups with %d obs per group:\n", n_groups, n_per_group))

      # Generate test data with shared and unique components
      set.seed(42)
      dataList <- generate_hierarchical_beta_data(n_groups, n_per_group)

      # Set up priors
      priorParameters <- c(2, 8)
      hyperPriorParameters <- c(1, 0.125)
      gammaPriors <- c(2, 4)
      alphaPriors <- c(2, 4)
      mhStepSize <- c(0.1, 0.1)
      numSticks <- 10

      # Benchmark R implementation
      times_r <- numeric(n_reps)
      for (i in 1:n_reps) {
        enable_cpp_hierarchical_samplers(FALSE)
        time_r <- system.time({
          dp_r <- DirichletProcessHierarchicalBeta(
            dataList, maxY = 1,
            priorParameters = priorParameters,
            hyperPriorParameters = hyperPriorParameters,
            gammaPriors = gammaPriors,
            alphaPriors = alphaPriors,
            mhStepSize = mhStepSize,
            numSticks = numSticks
          )
          dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
        })
        times_r[i] <- time_r["elapsed"]
      }

      # Benchmark C++ implementation
      times_cpp <- numeric(n_reps)
      global_clusters_cpp <- numeric(n_reps)

      if (exists("hierarchical_beta_fit_cpp")) {
        for (i in 1:n_reps) {
          enable_cpp_hierarchical_samplers(TRUE)
          time_cpp <- system.time({
            dp_cpp <- DirichletProcessHierarchicalBeta(
              dataList, maxY = 1,
              priorParameters = priorParameters,
              hyperPriorParameters = hyperPriorParameters,
              gammaPriors = gammaPriors,
              alphaPriors = alphaPriors,
              mhStepSize = mhStepSize,
              numSticks = numSticks
            )
            dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
          })
          times_cpp[i] <- time_cpp["elapsed"]

          # Count unique global parameters
          global_params <- unique(unlist(lapply(dp_cpp$indDP, function(x) x$clusterParameters)))
          global_clusters_cpp[i] <- length(global_params) / 2  # Divide by 2 for (mu, tau) pairs
        }
      } else {
        times_cpp <- rep(NA, n_reps)
        global_clusters_cpp <- rep(NA, n_reps)
      }

      # Calculate statistics
      mean_r <- mean(times_r)
      mean_cpp <- mean(times_cpp, na.rm = TRUE)
      speedup <- if(!is.na(mean_cpp)) mean_r / mean_cpp else NA

      # Print results
      cat(sprintf("  R: %.3fs, C++: %.3fs, Speedup: %.1fx, Global Clusters: %.1f\n",
                  mean_r, mean_cpp, speedup, mean(global_clusters_cpp, na.rm = TRUE)))

      # Store results
      results <- rbind(results, data.frame(
        n_groups = n_groups,
        n_per_group = n_per_group,
        implementation = "R",
        time = times_r,
        memory_mb = NA,
        global_clusters_found = NA
      ))

      if (!all(is.na(times_cpp))) {
        results <- rbind(results, data.frame(
          n_groups = n_groups,
          n_per_group = n_per_group,
          implementation = "C++",
          time = times_cpp,
          memory_mb = NA,
          global_clusters_found = global_clusters_cpp
        ))
      }
    }
  }

  return(results)
}

#' Profile Hierarchical Beta components
#'
#' @param n_groups Number of groups
#' @param n_per_group Observations per group
#' @export
profile_hierarchical_beta_components <- function(n_groups = 3, n_per_group = 200) {

  cat("\n=== Hierarchical Beta Component Profiling ===\n\n")

  # Ensure gtools is loaded and attached
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' is required but not installed.")
  }
  library(gtools)  # Attach to search path for C++ access

  # Generate test data
  set.seed(123)
  dataList <- generate_hierarchical_beta_data(n_groups, n_per_group)

  # Set up parameters
  priorParameters <- c(2, 8)
  hyperPriorParameters <- c(1, 0.125)
  gammaPriors <- c(2, 4)
  alphaPriors <- c(2, 4)
  mhStepSize <- c(0.1, 0.1)
  numSticks <- 10

  # Initialize hierarchical DP
  cat("Initializing hierarchical DP...\n")
  enable_cpp_hierarchical_samplers(FALSE)
  dp_r <- DirichletProcessHierarchicalBeta(
    dataList, maxY = 1,
    priorParameters = priorParameters,
    hyperPriorParameters = hyperPriorParameters,
    gammaPriors = gammaPriors,
    alphaPriors = alphaPriors,
    mhStepSize = mhStepSize,
    numSticks = numSticks
  )

  if (exists("hierarchical_beta_fit_cpp")) {
    enable_cpp_hierarchical_samplers(TRUE)
    dp_cpp <- DirichletProcessHierarchicalBeta(
      dataList, maxY = 1,
      priorParameters = priorParameters,
      hyperPriorParameters = hyperPriorParameters,
      gammaPriors = gammaPriors,
      alphaPriors = alphaPriors,
      mhStepSize = mhStepSize,
      numSticks = numSticks
    )
  } else {
    cat("C++ implementation not available, skipping C++ profiling\n")
    # Return R-only results
    time_global_r <- system.time({
      for (i in 1:50) {
        dp_temp <- GlobalParameterUpdate(dp_r)
      }
    })["elapsed"] / 50 * 1000

    time_g0_r <- system.time({
      for (i in 1:50) {
        dp_temp <- UpdateG0(dp_r)
      }
    })["elapsed"] / 50 * 1000

    time_gamma_r <- system.time({
      for (i in 1:50) {
        dp_temp <- UpdateGamma(dp_r)
      }
    })["elapsed"] / 50 * 1000

    return(data.frame(
      component = c("GlobalParameterUpdate", "UpdateG0", "UpdateGamma"),
      implementation = "R",
      mean_time_ms = c(time_global_r, time_g0_r, time_gamma_r)
    ))
  }

  # Time individual components
  n_timing_reps <- 50

  cat("Profiling R implementation...\n")

  # Global parameter update
  time_global_r <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- GlobalParameterUpdate(dp_r)
    }
  })["elapsed"] / n_timing_reps * 1000  # Convert to ms

  # Update G0
  time_g0_r <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- UpdateG0(dp_r)
    }
  })["elapsed"] / n_timing_reps * 1000

  # Update gamma
  time_gamma_r <- system.time({
    for (i in 1:n_timing_reps) {
      dp_temp <- UpdateGamma(dp_r)
    }
  })["elapsed"] / n_timing_reps * 1000

  if (exists("hierarchical_beta_fit_cpp")) {
    cat("Profiling C++ implementation...\n")

    # Use tryCatch to handle potential C++ errors gracefully
    time_global_cpp <- tryCatch({
      system.time({
        for (i in 1:n_timing_reps) {
          dp_temp <- hierarchical_beta_global_parameter_update_cpp(dp_cpp)
        }
      })["elapsed"] / n_timing_reps * 1000
    }, error = function(e) {
      cat("Error in C++ global parameter update:", e$message, "\n")
      NA
    })

    time_g0_cpp <- tryCatch({
      system.time({
        for (i in 1:n_timing_reps) {
          dp_temp <- hierarchical_beta_update_g0_cpp(dp_cpp)
        }
      })["elapsed"] / n_timing_reps * 1000
    }, error = function(e) {
      cat("Error in C++ G0 update:", e$message, "\n")
      NA
    })

    time_gamma_cpp <- tryCatch({
      system.time({
        for (i in 1:n_timing_reps) {
          dp_temp <- hierarchical_beta_update_gamma_cpp(dp_cpp)
        }
      })["elapsed"] / n_timing_reps * 1000
    }, error = function(e) {
      cat("Error in C++ gamma update:", e$message, "\n")
      NA
    })
  } else {
    time_global_cpp <- time_g0_cpp <- time_gamma_cpp <- NA
  }

  # Create results dataframe
  component_results <- data.frame(
    component = c("GlobalParameterUpdate", "UpdateG0", "UpdateGamma"),
    implementation = rep(c("R", "C++"), each = 3),
    mean_time_ms = c(time_global_r, time_g0_r, time_gamma_r,
                     time_global_cpp, time_g0_cpp, time_gamma_cpp)
  )

  cat("\nComponent timings (milliseconds):\n")
  print(component_results)

  # Calculate speedup only if C++ results are available
  if (!all(is.na(c(time_global_cpp, time_g0_cpp, time_gamma_cpp)))) {
    speedup_df <- component_results %>%
      pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
      mutate(speedup = R / `C++`)

    cat("\nSpeedup by component:\n")
    print(speedup_df)
  }

  return(component_results)
}

#' Test statistical equivalence between R and C++ implementations
#'
#' @param n_groups Number of groups
#' @param n_per_group Observations per group
#' @param n_iter Number of iterations
#' @export
test_hierarchical_beta_statistical_equivalence <- function(n_groups = 2,
                                                           n_per_group = 100,
                                                           n_iter = 200) {

  cat("\n=== Testing Statistical Equivalence ===\n\n")

  if (!exists("hierarchical_beta_fit_cpp")) {
    cat("C++ implementation not available for comparison.\n")
    return(NULL)
  }

  # Generate test data
  set.seed(456)
  dataList <- generate_hierarchical_beta_data(n_groups, n_per_group)

  # Common parameters
  params <- list(
    maxY = 1,
    priorParameters = c(2, 8),
    hyperPriorParameters = c(1, 0.125),
    gammaPriors = c(2, 4),
    alphaPriors = c(2, 4),
    mhStepSize = c(0.1, 0.1),
    numSticks = 10
  )

  # Fit with R implementation
  set.seed(789)
  enable_cpp_hierarchical_samplers(FALSE)
  dp_r <- do.call(DirichletProcessHierarchicalBeta, c(list(dataList), params))
  dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)

  # Fit with C++ implementation
  set.seed(789)
  enable_cpp_hierarchical_samplers(TRUE)
  dp_cpp <- do.call(DirichletProcessHierarchicalBeta, c(list(dataList), params))
  dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)

  # Compare results
  cat("Number of clusters per group:\n")
  for (i in 1:n_groups) {
    cat(sprintf("  Group %d - R: %d, C++: %d\n",
                i, dp_r$indDP[[i]]$numberClusters, dp_cpp$indDP[[i]]$numberClusters))
  }

  # Compare gamma values
  cat("\nGamma parameter:\n")
  cat(sprintf("  R: %.3f\n", dp_r$gamma))
  cat(sprintf("  C++: %.3f\n", dp_cpp$gamma))

  # Compare alpha values
  cat("\nAlpha parameters:\n")
  for (i in 1:n_groups) {
    cat(sprintf("  Group %d - R: %.3f, C++: %.3f\n",
                i, dp_r$indDP[[i]]$alpha, dp_cpp$indDP[[i]]$alpha))
  }

  # Statistical test for parameter equivalence
  all_params_r <- unlist(lapply(dp_r$indDP, function(x) x$clusterParameters))
  all_params_cpp <- unlist(lapply(dp_cpp$indDP, function(x) x$clusterParameters))

  if (length(all_params_r) == length(all_params_cpp)) {
    ks_test <- ks.test(all_params_r, all_params_cpp)
    cat(sprintf("\nKS test for parameter distributions p-value: %.3f\n", ks_test$p.value))
  }

  return(list(dp_r = dp_r, dp_cpp = dp_cpp))
}

#' Comprehensive Hierarchical Beta benchmarking
#'
#' @param n_groups_vec Vector of number of groups
#' @param n_per_group_vec Vector of observations per group
#' @param n_iter_vec Vector of iteration counts
#' @param n_shared_clusters_vec Vector of shared cluster counts
#' @param n_reps Number of repetitions
#' @export
benchmark_hierarchical_beta_comprehensive <- function(
    n_groups_vec = c(2, 3, 5),
    n_per_group_vec = c(100, 200, 500),
    n_iter_vec = c(100, 250),
    n_shared_clusters_vec = c(1, 2),
    n_reps = 3) {

  cat("\n=== Comprehensive Hierarchical Beta Benchmarking ===\n\n")

  total_scenarios <- length(n_groups_vec) * length(n_per_group_vec) *
    length(n_iter_vec) * length(n_shared_clusters_vec)
  cat(sprintf("Testing %d scenarios with %d reps each...\n", total_scenarios, n_reps))

  results <- data.frame()
  scenario_id <- 0

  for (n_groups in n_groups_vec) {
    for (n_per_group in n_per_group_vec) {
      for (n_iter in n_iter_vec) {
        for (n_shared in n_shared_clusters_vec) {

          scenario_id <- scenario_id + 1
          cat(sprintf("\nScenario %d/%d: groups=%d, n_per_group=%d, iter=%d, shared=%d\n",
                      scenario_id, total_scenarios, n_groups, n_per_group, n_iter, n_shared))

          # Generate data with specified structure
          set.seed(scenario_id)

          # Create shared parameters
          shared_params <- list()
          for (s in 1:n_shared) {
            mu <- runif(1, 0.2, 0.8)
            tau <- runif(1, 5, 20)
            shared_params[[s]] <- c(mu, tau)
          }

          # Create unique parameters for each group
          unique_params <- list()
          for (g in 1:n_groups) {
            mu <- runif(1, 0.1, 0.9)
            tau <- runif(1, 5, 20)
            unique_params[[g]] <- c(mu, tau)
          }

          # Generate hierarchical data
          dataList <- list()
          for (g in 1:n_groups) {
            # Mix shared and unique components
            n_shared_obs <- round(n_per_group * 0.6)  # 60% shared
            n_unique_obs <- n_per_group - n_shared_obs

            data_shared <- numeric(0)
            if (n_shared > 0) {
              # Sample from shared components
              for (s in 1:n_shared) {
                n_s <- round(n_shared_obs / n_shared)
                params <- shared_params[[s]]
                data_shared <- c(data_shared,
                                 rbeta(n_s, params[1] * params[2],
                                       (1 - params[1]) * params[2]))
              }
            }

            # Sample from unique component
            params <- unique_params[[g]]
            data_unique <- rbeta(n_unique_obs, params[1] * params[2],
                                 (1 - params[1]) * params[2])

            dataList[[g]] <- c(data_shared, data_unique)
          }

          # Common parameters
          params <- list(
            maxY = 1,
            priorParameters = c(2, 8),
            hyperPriorParameters = c(1, 0.125),
            gammaPriors = c(2, 4),
            alphaPriors = c(2, 4),
            mhStepSize = c(0.1, 0.1),
            numSticks = 20
          )

          # Run R implementation
          for (rep in 1:n_reps) {
            enable_cpp_hierarchical_samplers(FALSE)

            time_r <- system.time({
              dp_r <- do.call(DirichletProcessHierarchicalBeta,
                              c(list(dataList), params))
              dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
            })["elapsed"]

            # Count clusters
            clusters_per_group <- sapply(dp_r$indDP, function(x) x$numberClusters)

            results <- rbind(results, data.frame(
              scenario_id = scenario_id,
              n_groups = n_groups,
              n_per_group = n_per_group,
              n_iter = n_iter,
              n_shared_true = n_shared,
              implementation = "R",
              rep = rep,
              time = time_r,
              gamma = dp_r$gamma,
              mean_alpha = mean(sapply(dp_r$indDP, function(x) x$alpha)),
              mean_clusters_per_group = mean(clusters_per_group),
              total_clusters = sum(clusters_per_group),
              memory_mb = NA
            ))
          }

          # Run C++ implementation if available
          if (exists("hierarchical_beta_fit_cpp")) {
            for (rep in 1:n_reps) {
              enable_cpp_hierarchical_samplers(TRUE)

              # Clear memory tracking if available
              if (exists("clear_memory_tracking")) {
                clear_memory_tracking()
              }

              time_cpp <- system.time({
                dp_cpp <- do.call(DirichletProcessHierarchicalBeta,
                                  c(list(dataList), params))
                dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
              })["elapsed"]

              # Get memory usage if available
              memory_mb <- NA
              if (exists("get_memory_tracking")) {
                mem_track <- get_memory_tracking()
                if (nrow(mem_track) > 0) {
                  memory_mb <- sum(mem_track$mb)
                }
              }

              # Count clusters
              clusters_per_group <- sapply(dp_cpp$indDP, function(x) x$numberClusters)

              results <- rbind(results, data.frame(
                scenario_id = scenario_id,
                n_groups = n_groups,
                n_per_group = n_per_group,
                n_iter = n_iter,
                n_shared_true = n_shared,
                implementation = "C++",
                rep = rep,
                time = time_cpp,
                gamma = dp_cpp$gamma,
                mean_alpha = mean(sapply(dp_cpp$indDP, function(x) x$alpha)),
                mean_clusters_per_group = mean(clusters_per_group),
                total_clusters = sum(clusters_per_group),
                memory_mb = memory_mb
              ))
            }

            # Print speedup for this scenario
            current_results <- results %>%
              filter(scenario_id == !!scenario_id) %>%
              group_by(implementation) %>%
              summarise(mean_time = mean(time), .groups = "drop")

            if (nrow(current_results) == 2) {
              speedup <- current_results$mean_time[current_results$implementation == "R"] /
                current_results$mean_time[current_results$implementation == "C++"]

              cat(sprintf("  Mean times - R: %.3fs, C++: %.3fs, Speedup: %.1fx\n",
                          current_results$mean_time[current_results$implementation == "R"],
                          current_results$mean_time[current_results$implementation == "C++"],
                          speedup))
            }
          }
        }
      }
    }
  }

  return(results)
}

#' Create visualization dashboard for Hierarchical Beta benchmarks
#'
#' @param bench_results Results from benchmark_hierarchical_beta_comprehensive
#' @export
visualize_hierarchical_beta_benchmarks <- function(bench_results) {

  # Check if we have C++ results
  has_cpp <- "C++" %in% bench_results$implementation

  # Set color scheme
  impl_colors <- c("R" = "#E41A1C", "C++" = "#377EB8")

  # 1. Speedup by number of groups and data size
  if (has_cpp) {
    speedup_data <- bench_results %>%
      group_by(n_groups, n_per_group, implementation) %>%
      summarise(mean_time = mean(time), .groups = "drop") %>%
      pivot_wider(names_from = implementation, values_from = mean_time) %>%
      mutate(speedup = R / `C++`)

    p_speedup <- ggplot(speedup_data, aes(x = n_per_group, y = speedup,
                                          color = factor(n_groups))) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = "C++ Speedup by Data Size and Number of Groups",
        x = "Observations per Group",
        y = "Speedup Factor (R time / C++ time)",
        color = "Number of\nGroups"
      ) +
      theme_minimal() +
      geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5)
  } else {
    p_speedup <- ggplot() +
      ggtitle("C++ Implementation Not Available") +
      theme_minimal()
  }

  # 2. Scaling behavior
  scaling_data <- bench_results %>%
    filter(n_iter == min(n_iter)) %>%
    group_by(n_groups, n_per_group, implementation) %>%
    summarise(mean_time = mean(time), .groups = "drop")

  p_scaling <- ggplot(scaling_data, aes(x = n_per_group, y = mean_time,
                                        color = implementation,
                                        linetype = factor(n_groups))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Computation Time Scaling",
      x = "Observations per Group",
      y = "Time (seconds)",
      color = "Implementation",
      linetype = "Number of\nGroups"
    ) +
    theme_minimal()

  # 3. Group complexity impact
  group_impact <- bench_results %>%
    filter(n_per_group == median(unique(n_per_group))) %>%
    group_by(n_groups, implementation) %>%
    summarise(
      mean_time = mean(time),
      se_time = sd(time) / sqrt(n()),
      .groups = "drop"
    )

  p_groups <- ggplot(group_impact, aes(x = n_groups, y = mean_time,
                                       fill = implementation)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                  position = position_dodge(0.8), width = 0.25) +
    scale_fill_manual(values = impl_colors) +
    labs(
      title = "Impact of Number of Groups",
      x = "Number of Groups",
      y = "Time (seconds)",
      fill = "Implementation"
    ) +
    theme_minimal()

  # 4. Clustering accuracy (gamma parameter recovery)
  gamma_recovery <- bench_results %>%
    group_by(n_groups, n_shared_true, implementation) %>%
    summarise(
      mean_gamma = mean(gamma),
      sd_gamma = sd(gamma),
      .groups = "drop"
    )

  p_gamma <- ggplot(gamma_recovery,
                    aes(x = n_shared_true, y = mean_gamma,
                        color = implementation, shape = factor(n_groups))) +
    geom_point(size = 3, position = position_dodge(0.3)) +
    geom_errorbar(aes(ymin = mean_gamma - sd_gamma, ymax = mean_gamma + sd_gamma),
                  position = position_dodge(0.3), width = 0.2) +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Gamma Parameter Recovery",
      x = "True Number of Shared Clusters",
      y = "Mean Estimated Gamma",
      color = "Implementation",
      shape = "Number of\nGroups"
    ) +
    theme_minimal()

  # 5. Memory usage (if available)
  if (has_cpp && !all(is.na(bench_results$memory_mb))) {
    memory_data <- bench_results %>%
      filter(implementation == "C++", !is.na(memory_mb)) %>%
      group_by(n_groups, n_per_group) %>%
      summarise(
        mean_memory = mean(memory_mb),
        se_memory = sd(memory_mb) / sqrt(n()),
        .groups = "drop"
      )

    p_memory <- ggplot(memory_data, aes(x = n_per_group, y = mean_memory,
                                        color = factor(n_groups))) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = mean_memory - se_memory,
                        ymax = mean_memory + se_memory),
                    width = 0.1) +
      scale_x_log10() +
      labs(
        title = "C++ Memory Usage",
        x = "Observations per Group",
        y = "Memory (MB)",
        color = "Number of\nGroups"
      ) +
      theme_minimal()
  } else {
    p_memory <- ggplot() +
      ggtitle("Memory Tracking Not Available") +
      theme_minimal()
  }

  # 6. Iterations impact
  iter_impact <- bench_results %>%
    filter(n_groups == median(unique(n_groups)),
           n_per_group == median(unique(n_per_group))) %>%
    group_by(n_iter, implementation) %>%
    summarise(
      mean_time = mean(time),
      .groups = "drop"
    )

  p_iterations <- ggplot(iter_impact, aes(x = n_iter, y = mean_time,
                                          color = implementation)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = impl_colors) +
    labs(
      title = "Impact of MCMC Iterations",
      x = "Number of Iterations",
      y = "Time (seconds)",
      color = "Implementation"
    ) +
    theme_minimal()

  # Combine plots
  dashboard <- (p_speedup | p_scaling) /
    (p_groups | p_gamma) /
    (p_memory | p_iterations) +
    plot_annotation(
      title = "Hierarchical Beta Dirichlet Process: R vs C++ Performance",
      subtitle = sprintf("Based on %d benchmark scenarios",
                         length(unique(bench_results$scenario_id)))
    )

  return(list(
    dashboard = dashboard,
    speedup = p_speedup,
    scaling = p_scaling,
    groups = p_groups,
    gamma = p_gamma,
    memory = p_memory,
    iterations = p_iterations
  ))
}

#' Run complete Hierarchical Beta benchmark suite and generate report
#'
#' @export
run_hierarchical_beta_benchmark_report <- function() {

  cat("\n")
  cat("================================================\n")
  cat("Hierarchical Beta Dirichlet Process Benchmarking\n")
  cat("================================================\n")

  # Ensure required packages are loaded
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' is required but not installed.")
  }
  library(gtools)  # Attach to search path

  # Suppress specific warnings about C++ fallback
  suppressWarnings({
    # Check if C++ is available
    cpp_available <- exists("hierarchical_beta_fit_cpp")
  })

  # Check if C++ is available
  cpp_available <- exists("hierarchical_beta_fit_cpp")

  if (!cpp_available) {
    cat("\n⚠️  WARNING: Hierarchical Beta C++ implementation not available.\n")
    cat("   Running R-only benchmarks.\n\n")
  } else {
    cat("\n✓ C++ implementation detected. Running comparative benchmarks.\n\n")
  }

  # 1. Quick benchmark
  cat("\n1. Running quick benchmark...\n")
  quick_results <- quick_hierarchical_beta_benchmark()

  # 2. Component profiling
  cat("\n2. Profiling individual components...\n")
  component_results <- profile_hierarchical_beta_components()

  # 3. Statistical equivalence
  if (cpp_available) {
    cat("\n3. Testing statistical equivalence...\n")
    equiv_results <- test_hierarchical_beta_statistical_equivalence()
  } else {
    cat("\n3. Skipping statistical equivalence test (C++ not available).\n")
    equiv_results <- NULL
  }

  # 4. Comprehensive benchmarks
  cat("\n4. Running comprehensive benchmarks...\n")
  bench_results <- benchmark_hierarchical_beta_comprehensive(
    n_groups_vec = c(2, 3, 5),
    n_per_group_vec = c(100, 200),
    n_iter_vec = c(100),
    n_shared_clusters_vec = c(1, 2),
    n_reps = 3
  )

  # Generate summary statistics
  summary_stats <- bench_results %>%
    group_by(implementation) %>%
    summarise(
      mean_time = mean(time),
      median_time = median(time),
      mean_memory = mean(memory_mb, na.rm = TRUE),
      mean_gamma = mean(gamma),
      .groups = "drop"
    )

  # Calculate speedup by scenario
  if (cpp_available) {
    speedup_by_scenario <- bench_results %>%
      group_by(n_groups, n_per_group, implementation) %>%
      summarise(mean_time = mean(time), .groups = "drop") %>%
      pivot_wider(names_from = implementation, values_from = mean_time) %>%
      mutate(speedup = R / `C++`)
  }

  # Print summary report
  cat("\n")
  cat("================================================\n")
  cat("HIERARCHICAL BETA BENCHMARK SUMMARY\n")
  cat("================================================\n")

  cat("\n1. OVERALL PERFORMANCE:\n")
  print(summary_stats)

  if (cpp_available) {
    cat("\n2. SPEEDUP BY SCENARIO:\n")
    print(speedup_by_scenario)

    cat("\n3. COMPONENT-LEVEL SPEEDUP:\n")
    comp_speedup <- component_results %>%
      pivot_wider(names_from = implementation, values_from = mean_time_ms) %>%
      mutate(speedup = R / `C++`)
    print(comp_speedup)

    cat("\n4. OVERALL SPEEDUP RANGE:\n")
    overall_speedup <- summary_stats$mean_time[summary_stats$implementation == "R"] /
      summary_stats$mean_time[summary_stats$implementation == "C++"]
    cat(sprintf("  Mean speedup: %.1fx\n", overall_speedup))
    cat(sprintf("  Min speedup: %.1fx\n", min(speedup_by_scenario$speedup)))
    cat(sprintf("  Max speedup: %.1fx\n", max(speedup_by_scenario$speedup)))
  }

  # Create visualizations
  cat("\n5. Generating visualizations...\n")
  plots <- visualize_hierarchical_beta_benchmarks(bench_results)

  # Display main dashboard
  print(plots$dashboard)

  # Save results
  results <- list(
    quick_benchmark = quick_results,
    component_profile = component_results,
    statistical_equivalence = equiv_results,
    comprehensive_benchmark = bench_results,
    summary_stats = summary_stats,
    speedup_by_scenario = if(cpp_available) speedup_by_scenario else NULL,
    plots = plots
  )

  cat("\nBenchmarking complete!\n")

  return(invisible(results))
}

# Run the benchmark if sourced interactively
if (interactive()) {
  results <- run_hierarchical_beta_benchmark_report()
}
