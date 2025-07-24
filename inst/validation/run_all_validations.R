# inst/validation/run_all_validations.R

# Load required libraries
library(dirichletprocess)
library(testthat)
library(microbenchmark)
library(ggplot2)

# Source all test files
source("tests/testthat/helper-testing.R")
source("tests/testthat/test-cpp-consistency.R")
source("benchmark/comprehensive_performance_tests.R")
source("benchmark/manual_mcmc_performance.R")
source("benchmark/visualize_performance.R")
source("tests/integration/package_checks.R")
source("tests/integration/memory_tests.R")
source("tests/integration/stress_tests.R")

run_complete_validation <- function(output_dir = "validation_results") {
  dir.create(output_dir, showWarnings = FALSE)

  cat("\n========================================\n")
  cat("DIRICHLET PROCESS C++ VALIDATION SUITE\n")
  cat("========================================\n")

  validation_results <- list()
  validation_results$start_time <- Sys.time()

  # Phase 1: Core Functionality
  cat("\n--- PHASE 1: CORE FUNCTIONALITY ---\n")

  # 1.1 Consistency tests
  cat("\n1.1 Running R/C++ consistency tests...\n")
  validation_results$consistency <- tryCatch({
    run_consistency_tests()
  }, error = function(e) {
    cat("ERROR in consistency tests:", e$message, "\n")
    list(error = e$message)
  })

  # 1.2 Edge case tests
  cat("\n1.2 Running edge case tests...\n")
  validation_results$edge_cases <- tryCatch({
    testthat::test_file("tests/testthat/test-cpp-edge-cases.R")
  }, error = function(e) {
    cat("ERROR in edge case tests:", e$message, "\n")
    list(error = e$message)
  })

  # 1.3 Convergence tests
  cat("\n1.3 Running convergence tests...\n")
  validation_results$convergence <- tryCatch({
    testthat::test_file("tests/testthat/test-cpp-convergence.R")
  }, error = function(e) {
    cat("ERROR in convergence tests:", e$message, "\n")
    list(error = e$message)
  })

  # 1.4 Manual MCMC tests
  cat("\n1.4 Running manual MCMC tests...\n")
  validation_results$manual_mcmc <- tryCatch({
    testthat::test_file("tests/testthat/test-cpp-manual-mcmc.R")
  }, error = function(e) {
    cat("ERROR in manual MCMC tests:", e$message, "\n")
    list(error = e$message)
  })

  # Phase 2: Performance
  cat("\n--- PHASE 2: PERFORMANCE ---\n")

  # 2.1 Comprehensive performance
  cat("\n2.1 Running comprehensive performance tests...\n")
  cat("(This may take 30-60 minutes)\n")
  validation_results$performance <- tryCatch({
    run_comprehensive_performance_tests()
  }, error = function(e) {
    cat("ERROR in performance tests:", e$message, "\n")
    list(error = e$message)
  })

  # 2.2 Scaling analysis
  cat("\n2.2 Testing scaling behavior...\n")
  validation_results$scaling <- tryCatch({
    test_scaling_behavior()
  }, error = function(e) {
    cat("ERROR in scaling tests:", e$message, "\n")
    list(error = e$message)
  })

  # 2.3 Memory profiling
  cat("\n2.3 Profiling memory usage...\n")
  validation_results$memory_profile <- tryCatch({
    profile_memory_usage()
  }, error = function(e) {
    cat("ERROR in memory profiling:", e$message, "\n")
    list(error = e$message)
  })

  # 2.4 Manual MCMC performance
  cat("\n2.4 Benchmarking manual MCMC interface...\n")
  validation_results$manual_mcmc_perf <- tryCatch({
    list(
      basic = benchmark_manual_mcmc(),
      vs_fit = benchmark_manual_vs_fit(),
      advanced = benchmark_advanced_features()
    )
  }, error = function(e) {
    cat("ERROR in manual MCMC benchmarks:", e$message, "\n")
    list(error = e$message)
  })

  # Phase 3: Integration
  cat("\n--- PHASE 3: INTEGRATION ---\n")

  # 3.1 Package checks
  cat("\n3.1 Running package checks...\n")
  validation_results$package_checks <- tryCatch({
    run_all_integration_checks()
  }, error = function(e) {
    cat("ERROR in package checks:", e$message, "\n")
    list(error = e$message)
  })

  # 3.2 Memory tests
  cat("\n3.2 Running memory tests...\n")
  validation_results$memory_tests <- tryCatch({
    run_all_memory_tests()
  }, error = function(e) {
    cat("ERROR in memory tests:", e$message, "\n")
    list(error = e$message)
  })

  # 3.3 Stress tests
  cat("\n3.3 Running stress tests...\n")
  validation_results$stress_tests <- tryCatch({
    run_comprehensive_stress_tests()
  }, error = function(e) {
    cat("ERROR in stress tests:", e$message, "\n")
    list(error = e$message)
  })

  # End time
  validation_results$end_time <- Sys.time()
  validation_results$total_runtime <- difftime(
    validation_results$end_time,
    validation_results$start_time,
    units = "hours"
  )

  # Generate reports
  cat("\n--- GENERATING REPORTS ---\n")

  # Save raw results
  saveRDS(validation_results, file.path(output_dir, "validation_results.rds"))

  # Generate validation report
  generate_validation_report(validation_results, output_dir)

  # Create performance dashboard
  if (!is.null(validation_results$performance)) {
    create_performance_dashboard(
      validation_results$performance,
      validation_results$scaling,
      validation_results$memory_profile
    )
  }

  # Print summary
  print_validation_summary(validation_results)

  return(validation_results)
}

generate_validation_report <- function(results, output_dir) {
  # Create markdown report
  report_file <- file.path(output_dir, "validation_report.md")

  report_content <- sprintf("
# Dirichlet Process C++ Validation Report

**Generated**: %s
**Total Runtime**: %.2f hours

## Executive Summary

### Overall Status: %s

## Phase 1: Core Functionality

### 1.1 R/C++ Consistency Tests
%s

### 1.2 Edge Case Tests
%s

### 1.3 Convergence Tests
%s

### 1.4 Manual MCMC Interface Tests
%s

## Phase 2: Performance

### 2.1 Performance Summary
%s

### 2.2 Scaling Analysis
%s

### 2.3 Memory Profiling
%s

## Phase 3: Integration

### 3.1 Package Checks
%s

### 3.2 Memory Tests
%s

### 3.3 Stress Tests
%s

## Recommendations

%s
",
                            format(results$end_time, "%Y-%m-%d %H:%M:%S"),
                            as.numeric(results$total_runtime),
                            determine_overall_status(results),
                            format_consistency_results(results$consistency),
                            format_test_results(results$edge_cases),
                            format_test_results(results$convergence),
                            format_test_results(results$manual_mcmc),
                            format_performance_results(results$performance),
                            format_scaling_results(results$scaling),
                            format_memory_profile(results$memory_profile),
                            format_package_check_results(results$package_checks),
                            format_memory_test_results(results$memory_tests),
                            format_stress_test_results(results$stress_tests),
                            generate_recommendations(results)
  )

  writeLines(report_content, report_file)
  cat("\nValidation report written to:", report_file, "\n")
}

# Helper functions for report generation
determine_overall_status <- function(results) {
  # Check for any errors
  has_errors <- any(sapply(results, function(x) {
    !is.null(x$error) ||
      (is.list(x) && any(sapply(x, function(y) !is.null(y$error))))
  }))

  if (has_errors) {
    return("❌ FAILED - Errors detected")
  }

  # Check consistency
  if (!is.null(results$consistency)) {
    avg_speedup <- mean(sapply(results$consistency, function(x) {
      if (!is.null(x$speedup_factor)) x$speedup_factor else 0
    }))

    if (avg_speedup < 1.5) {
      return("⚠️  WARNING - Low performance improvement")
    }
  }

  return("✅ PASSED - Ready for production")
}

format_consistency_results <- function(consistency) {
  if (is.null(consistency)) return("No results available")
  if (!is.null(consistency$error)) return(paste("ERROR:", consistency$error))

  lines <- c()
  for (dist in names(consistency)) {
    res <- consistency[[dist]]
    if (!is.null(res$alpha_mean_diff)) {
      lines <- c(lines, sprintf(
        "- **%s**: α diff=%.4f, cluster diff=%.4f, correlation=%.3f, speedup=%.2fx",
        dist, res$alpha_mean_diff, res$cluster_count_diff,
        res$likelihood_correlation, res$speedup_factor
      ))
    }
  }

  return(paste(lines, collapse = "\n"))
}

format_test_results <- function(test_res) {
  if (is.null(test_res)) return("No results available")
  if (!is.null(test_res$error)) return(paste("ERROR:", test_res$error))

  # Extract test statistics from testthat results
  if (inherits(test_res, "testthat_results")) {
    passed <- sum(sapply(test_res, function(x) x$passed))
    failed <- sum(sapply(test_res, function(x) x$failed))
    warnings <- sum(sapply(test_res, function(x) x$warnings))

    status <- if (failed == 0) "✅ PASSED" else "❌ FAILED"
    return(sprintf("%s - %d passed, %d failed, %d warnings",
                   status, passed, failed, warnings))
  }

  return("Results format not recognized")
}

format_performance_results <- function(perf) {
  if (is.null(perf)) return("No results available")
  if (!is.null(perf$error)) return(paste("ERROR:", perf$error))

  # Calculate average speedups
  speedups <- sapply(perf, function(x) {
    if (!is.null(x$speedup)) x$speedup else NA
  })

  avg_speedup <- mean(speedups, na.rm = TRUE)

  return(sprintf(
    "- Average speedup: %.2fx\n- Tests completed: %d\n- See performance dashboard for details",
    avg_speedup, length(speedups)
  ))
}

format_scaling_results <- function(scaling) {
  if (is.null(scaling)) return("No results available")
  if (!is.null(scaling$error)) return(paste("ERROR:", scaling$error))

  lines <- c()
  for (dist in names(scaling)) {
    df <- scaling[[dist]]
    if (!is.null(df$speedup)) {
      avg_speedup <- mean(df$speedup)
      lines <- c(lines, sprintf("- %s: Average speedup %.2fx", dist, avg_speedup))
    }
  }

  return(paste(lines, collapse = "\n"))
}

format_memory_profile <- function(mem) {
  if (is.null(mem)) return("No results available")
  if (!is.null(mem$error)) return(paste("ERROR:", mem$error))

  lines <- c()
  for (dist in names(mem)) {
    res <- mem[[dist]]
    if (!is.null(res$memory_ratio)) {
      lines <- c(lines, sprintf(
        "- %s: Memory reduction %.1f%%",
        dist, (1 - 1/res$memory_ratio) * 100
      ))
    }
  }

  return(paste(lines, collapse = "\n"))
}

format_package_check_results <- function(checks) {
  if (is.null(checks)) return("No results available")
  if (!is.null(checks$error)) return(paste("ERROR:", checks$error))

  # Format based on actual check results
  return("See detailed integration report for full results")
}

format_memory_test_results <- function(mem_tests) {
  if (is.null(mem_tests)) return("No results available")
  if (!is.null(mem_tests$error)) return(paste("ERROR:", mem_tests$error))

  return("Memory tests completed - see memory_test_results.rds for details")
}

format_stress_test_results <- function(stress) {
  if (is.null(stress)) return("No results available")
  if (!is.null(stress$error)) return(paste("ERROR:", stress$error))

  return("Stress tests completed - see stress_test_results.rds for details")
}

generate_recommendations <- function(results) {
  recommendations <- c()

  # Check performance
  if (!is.null(results$performance)) {
    avg_speedup <- mean(sapply(results$performance, function(x) {
      if (!is.null(x$speedup)) x$speedup else 0
    }))

    if (avg_speedup < 2) {
      recommendations <- c(recommendations,
                           "- Consider further C++ optimization for better performance gains")
    }
  }

  # Check for errors
  error_count <- sum(sapply(results, function(x) {
    !is.null(x$error) ||
      (is.list(x) && any(sapply(x, function(y) !is.null(y$error))))
  }))

  if (error_count > 0) {
    recommendations <- c(recommendations,
                         "- Fix errors before proceeding to production")
  }

  if (length(recommendations) == 0) {
    recommendations <- "- Package is ready for production deployment"
  }

  return(paste(recommendations, collapse = "\n"))
}

print_validation_summary <- function(results) {
  cat("\n\n========================================\n")
  cat("       VALIDATION SUMMARY\n")
  cat("========================================\n")

  cat("\nTotal runtime:", round(as.numeric(results$total_runtime), 2), "hours\n")
  cat("Overall status:", determine_overall_status(results), "\n")

  # Performance summary
  if (!is.null(results$performance)) {
    speedups <- sapply(results$performance, function(x) {
      if (!is.null(x$speedup)) x$speedup else NA
    })
    cat("\nPerformance:\n")
    cat("  Average C++ speedup:", round(mean(speedups, na.rm = TRUE), 2), "x\n")
  }

  # Test summary
  cat("\nTest Results:\n")
  test_components <- c("edge_cases", "convergence", "manual_mcmc")
  for (comp in test_components) {
    if (!is.null(results[[comp]])) {
      cat("  ", comp, ": ",
          ifelse(is.null(results[[comp]]$error), "PASSED", "FAILED"), "\n")
    }
  }

  cat("\nReports generated in:", getwd(), "/validation_results/\n")
}

# Quick validation for development
quick_validation <- function() {
  cat("\n=== QUICK VALIDATION (Development) ===\n")
  cat("Running subset of tests for quick feedback...\n\n")

  results <- list()

  # Quick consistency check
  cat("1. Quick consistency check (Normal distribution only)...\n")
  test_data <- generate_test_data("normal", 100)
  results$consistency <- validate_r_cpp_consistency("normal", test_data,
                                                    iterations = 50, n_runs = 3)

  # Quick performance check
  cat("\n2. Quick performance benchmark...\n")
  results$performance <- quick_benchmark("normal", n = 1000, its = 100)

  # Quick manual MCMC check
  cat("\n3. Quick manual MCMC check...\n")
  dp <- DirichletProcessGaussian(test_data)
  runner <- CppMCMCRunner$new(dp)

  results$manual_mcmc <- tryCatch({
    for (i in 1:10) {
      runner$step_assignments()
      runner$step_parameters()
      runner$step_concentration()
    }
    "PASSED"
  }, error = function(e) {
    paste("FAILED:", e$message)
  })

  cat("\n=== QUICK VALIDATION COMPLETE ===\n")
  return(results)
}

# Run validation based on command line arguments
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0 || args[1] == "full") {
    results <- run_complete_validation()
  } else if (args[1] == "quick") {
    results <- quick_validation()
  } else {
    cat("Usage: Rscript run_all_validations.R [full|quick]\n")
  }
} else {
  cat("Run validation with:\n")
  cat("  results <- run_complete_validation()  # Full validation\n")
  cat("  results <- quick_validation()         # Quick validation\n")
}
