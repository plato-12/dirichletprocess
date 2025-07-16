# Practical MVNormal Benchmark - 20-30 Minutes Completion
# =====================================================
#
# Based on actual performance measurements showing 3-6 seconds per MCMC iteration,
# this benchmark is designed to complete in 20-30 minutes while providing
# meaningful performance comparisons across all covariance models.

library(dirichletprocess)
devtools::load_all()
set_use_cpp(TRUE)

cat("=== PRACTICAL MVNORMAL BENCHMARK ===\n")
cat("Designed for 20-30 minute completion with meaningful results\n\n")

# ==========================================
# OPTIMIZED CONFIGURATION
# ==========================================

# Configuration based on actual performance (5 seconds per MCMC iteration)
PRACTICAL_CONFIG <- list(
  mcmc_iterations = 5,      # 5 iterations = 25 seconds per test
  repetitions = 3,          # 3 repetitions for statistical reliability
  sample_sizes = c(50, 200), # 2 sample sizes for scalability testing
  dimensions = c(1, 2, 5),  # 3 dimensions for comprehensive coverage
  warmup_runs = 0           # No warmup to save time
)

# Calculate realistic timing
models_1d <- c("E", "V", "FULL")
models_2d <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")

combinations_1d <- length(models_1d) * length(PRACTICAL_CONFIG$sample_sizes)  # 3 * 2 = 6
combinations_2d <- length(models_2d) * length(PRACTICAL_CONFIG$sample_sizes) * 2  # 7 * 2 * 2 = 28 (for 2D and 5D)
total_combinations <- combinations_1d + combinations_2d  # 6 + 28 = 34

total_tests <- total_combinations * PRACTICAL_CONFIG$repetitions  # 34 * 3 = 102
estimated_time_minutes <- (total_tests * PRACTICAL_CONFIG$mcmc_iterations * 5) / 60  # 102 * 5 * 5 / 60 = 42.5 minutes

cat("Realistic estimates:\n")
cat("  Total combinations:", total_combinations, "\n")
cat("  Total tests:", total_tests, "\n")
cat("  Estimated time:", round(estimated_time_minutes, 1), "minutes\n")
cat("  Time per test: ~25 seconds\n\n")

# ==========================================
# BENCHMARK EXECUTION
# ==========================================

# Results storage
all_results <- list()
start_time <- Sys.time()
test_count <- 0

# Helper function to run single test
run_single_test <- function(model, dimension, sample_size, repetition) {
  tryCatch({
    # Generate appropriate data
    if (dimension == 1) {
      test_data <- matrix(rnorm(sample_size), ncol = 1)
    } else {
      # Create clustered data for better testing
      cluster_size <- sample_size %/% 3
      data1 <- mvtnorm::rmvnorm(cluster_size, mean = rep(-2, dimension), sigma = diag(dimension) * 0.5)
      data2 <- mvtnorm::rmvnorm(cluster_size, mean = rep(0, dimension), sigma = diag(dimension) * 0.5)
      data3 <- mvtnorm::rmvnorm(sample_size - 2*cluster_size, mean = rep(2, dimension), sigma = diag(dimension) * 0.5)
      test_data <- rbind(data1, data2, data3)
    }
    
    # Create prior parameters
    prior_params <- list(
      mu0 = rep(0, dimension),
      kappa0 = 1,
      nu = dimension + 2,
      Lambda = diag(dimension),
      covModel = model
    )
    
    # Run benchmark
    iter_start <- Sys.time()
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(test_data, md)
    dp <- Initialise(dp, numInitialClusters = 1)
    result <- Fit(dp, PRACTICAL_CONFIG$mcmc_iterations, progressBar = FALSE)
    iter_end <- Sys.time()
    
    # Return results
    return(list(
      success = TRUE,
      time_sec = as.numeric(iter_end - iter_start),
      time_per_mcmc = as.numeric(iter_end - iter_start) / PRACTICAL_CONFIG$mcmc_iterations,
      clusters = result$numberClusters,
      model = model,
      dimension = dimension,
      sample_size = sample_size,
      repetition = repetition
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      error = e$message,
      model = model,
      dimension = dimension,
      sample_size = sample_size,
      repetition = repetition
    ))
  })
}

# Test 1D models
cat("=== 1D MODELS ===\n")
for (model in models_1d) {
  for (n in PRACTICAL_CONFIG$sample_sizes) {
    cat(sprintf("Testing %s (d=1, n=%d):", model, n))
    
    model_results <- list()
    success_count <- 0
    
    for (rep in 1:PRACTICAL_CONFIG$repetitions) {
      test_count <- test_count + 1
      result <- run_single_test(model, 1, n, rep)
      
      if (result$success) {
        model_results[[rep]] <- result
        success_count <- success_count + 1
      } else {
        cat(sprintf(" [%d:ERROR]", rep))
      }
    }
    
    if (success_count > 0) {
      avg_time <- mean(sapply(model_results, function(x) x$time_sec))
      avg_clusters <- mean(sapply(model_results, function(x) x$clusters))
      avg_time_per_mcmc <- mean(sapply(model_results, function(x) x$time_per_mcmc))
      
      cat(sprintf(" ✓ %d/%d success, %.1fs total, %.1fs/iter, %d clusters\n", 
                  success_count, PRACTICAL_CONFIG$repetitions, avg_time, avg_time_per_mcmc, round(avg_clusters)))
      
      all_results <- append(all_results, model_results)
    } else {
      cat(" ✗ ALL FAILED\n")
    }
    
    # Progress indicator
    cat(sprintf("    Progress: %d/%d tests completed (%.1f%%)\n", test_count, total_tests, 100 * test_count / total_tests))
  }
}

# Test 2D models
cat("\n=== 2D MODELS ===\n")
for (model in models_2d) {
  for (n in PRACTICAL_CONFIG$sample_sizes) {
    cat(sprintf("Testing %s (d=2, n=%d):", model, n))
    
    model_results <- list()
    success_count <- 0
    
    for (rep in 1:PRACTICAL_CONFIG$repetitions) {
      test_count <- test_count + 1
      result <- run_single_test(model, 2, n, rep)
      
      if (result$success) {
        model_results[[rep]] <- result
        success_count <- success_count + 1
      } else {
        cat(sprintf(" [%d:ERROR]", rep))
      }
    }
    
    if (success_count > 0) {
      avg_time <- mean(sapply(model_results, function(x) x$time_sec))
      avg_clusters <- mean(sapply(model_results, function(x) x$clusters))
      avg_time_per_mcmc <- mean(sapply(model_results, function(x) x$time_per_mcmc))
      
      cat(sprintf(" ✓ %d/%d success, %.1fs total, %.1fs/iter, %d clusters\n", 
                  success_count, PRACTICAL_CONFIG$repetitions, avg_time, avg_time_per_mcmc, round(avg_clusters)))
      
      all_results <- append(all_results, model_results)
    } else {
      cat(" ✗ ALL FAILED\n")
    }
    
    # Progress indicator
    cat(sprintf("    Progress: %d/%d tests completed (%.1f%%)\n", test_count, total_tests, 100 * test_count / total_tests))
  }
}

# Test 5D models
cat("\n=== 5D MODELS ===\n")
for (model in models_2d) {  # Same models as 2D
  for (n in PRACTICAL_CONFIG$sample_sizes) {
    cat(sprintf("Testing %s (d=5, n=%d):", model, n))
    
    model_results <- list()
    success_count <- 0
    
    for (rep in 1:PRACTICAL_CONFIG$repetitions) {
      test_count <- test_count + 1
      result <- run_single_test(model, 5, n, rep)
      
      if (result$success) {
        model_results[[rep]] <- result
        success_count <- success_count + 1
      } else {
        cat(sprintf(" [%d:ERROR]", rep))
      }
    }
    
    if (success_count > 0) {
      avg_time <- mean(sapply(model_results, function(x) x$time_sec))
      avg_clusters <- mean(sapply(model_results, function(x) x$clusters))
      avg_time_per_mcmc <- mean(sapply(model_results, function(x) x$time_per_mcmc))
      
      cat(sprintf(" ✓ %d/%d success, %.1fs total, %.1fs/iter, %d clusters\n", 
                  success_count, PRACTICAL_CONFIG$repetitions, avg_time, avg_time_per_mcmc, round(avg_clusters)))
      
      all_results <- append(all_results, model_results)
    } else {
      cat(" ✗ ALL FAILED\n")
    }
    
    # Progress indicator
    cat(sprintf("    Progress: %d/%d tests completed (%.1f%%)\n", test_count, total_tests, 100 * test_count / total_tests))
  }
}

# ==========================================
# RESULTS ANALYSIS
# ==========================================

end_time <- Sys.time()
total_time <- end_time - start_time

cat("\n=== PRACTICAL BENCHMARK RESULTS ===\n")
cat("Total execution time:", format(total_time), "\n")
cat("Tests completed:", length(all_results), "\n")
cat("Success rate:", sprintf("%.1f%%", 100 * length(all_results) / total_tests), "\n")

if (length(all_results) > 0) {
  # Create summary dataframe
  summary_df <- data.frame(
    model = sapply(all_results, function(x) x$model),
    dimension = sapply(all_results, function(x) x$dimension),
    sample_size = sapply(all_results, function(x) x$sample_size),
    time_sec = sapply(all_results, function(x) x$time_sec),
    time_per_mcmc = sapply(all_results, function(x) x$time_per_mcmc),
    clusters = sapply(all_results, function(x) x$clusters),
    stringsAsFactors = FALSE
  )
  
  # Load dplyr for analysis
  library(dplyr)
  
  # Calculate model performance summary
  performance_summary <- summary_df %>%
    group_by(model, dimension, sample_size) %>%
    summarise(
      n_tests = n(),
      mean_time = mean(time_sec),
      sd_time = sd(time_sec),
      mean_time_per_mcmc = mean(time_per_mcmc),
      mean_clusters = mean(clusters),
      .groups = "drop"
    ) %>%
    arrange(mean_time_per_mcmc)
  
  cat("\n=== TOP 10 FASTEST MODELS ===\n")
  print(head(performance_summary, 10))
  
  cat("\n=== TOP 5 SLOWEST MODELS ===\n")
  print(tail(performance_summary, 5))
  
  # Performance by dimension
  cat("\n=== PERFORMANCE BY DIMENSION ===\n")
  dimension_summary <- summary_df %>%
    group_by(dimension) %>%
    summarise(
      models_tested = n_distinct(model),
      mean_time_per_mcmc = mean(time_per_mcmc),
      fastest_model = model[which.min(time_per_mcmc)],
      slowest_model = model[which.max(time_per_mcmc)],
      .groups = "drop"
    )
  
  print(dimension_summary)
  
  # Overall statistics
  cat("\n=== OVERALL STATISTICS ===\n")
  cat("Average time per MCMC iteration:", sprintf("%.2f seconds", mean(summary_df$time_per_mcmc)), "\n")
  cat("Fastest MCMC time:", sprintf("%.2f seconds", min(summary_df$time_per_mcmc)), "\n")
  cat("Slowest MCMC time:", sprintf("%.2f seconds", max(summary_df$time_per_mcmc)), "\n")
  cat("Average clusters found:", sprintf("%.1f", mean(summary_df$clusters)), "\n")
  
  # Model recommendations
  cat("\n=== MODEL RECOMMENDATIONS ===\n")
  
  # Best models by dimension
  best_1d <- performance_summary %>% filter(dimension == 1) %>% slice_min(mean_time_per_mcmc, n = 3)
  best_2d <- performance_summary %>% filter(dimension == 2) %>% slice_min(mean_time_per_mcmc, n = 3)
  best_5d <- performance_summary %>% filter(dimension == 5) %>% slice_min(mean_time_per_mcmc, n = 3)
  
  cat("Best 1D models:\n")
  for (i in 1:nrow(best_1d)) {
    cat(sprintf("  %d. %s: %.2f sec/iter\n", i, best_1d$model[i], best_1d$mean_time_per_mcmc[i]))
  }
  
  cat("Best 2D models:\n")
  for (i in 1:nrow(best_2d)) {
    cat(sprintf("  %d. %s: %.2f sec/iter\n", i, best_2d$model[i], best_2d$mean_time_per_mcmc[i]))
  }
  
  cat("Best 5D models:\n")
  for (i in 1:nrow(best_5d)) {
    cat(sprintf("  %d. %s: %.2f sec/iter\n", i, best_5d$model[i], best_5d$mean_time_per_mcmc[i]))
  }
  
  # Save results
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  save(all_results, summary_df, performance_summary, dimension_summary, 
       file = paste0("practical_benchmark_results_", timestamp, ".RData"))
  
  # Export CSV for external analysis
  write.csv(summary_df, paste0("practical_benchmark_raw_", timestamp, ".csv"), row.names = FALSE)
  write.csv(performance_summary, paste0("practical_benchmark_summary_", timestamp, ".csv"), row.names = FALSE)
  
  cat("\nResults saved to:\n")
  cat("  practical_benchmark_results_", timestamp, ".RData\n", sep = "")
  cat("  practical_benchmark_raw_", timestamp, ".csv\n", sep = "")
  cat("  practical_benchmark_summary_", timestamp, ".csv\n", sep = "")
  
} else {
  cat("No successful tests completed.\n")
}

cat("\n=== BENCHMARK COMPLETE ===\n")
cat("This benchmark provides practical performance comparisons\n")
cat("across all MVNormal covariance models in a reasonable timeframe.\n")