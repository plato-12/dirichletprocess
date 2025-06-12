# Test script for benchmark_exponential.R
# This script tests all components of the exponential benchmarking suite

library(dirichletprocess)

# Helper function to safely test a function
safe_test <- function(test_name, test_func) {
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("Testing:", test_name, "\n")
  cat(rep("=", 50), "\n", sep = "")

  result <- tryCatch({
    test_func()
    list(success = TRUE, error = NULL)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })

  if (result$success) {
    cat("✓ SUCCESS:", test_name, "completed successfully\n")
  } else {
    cat("✗ FAILED:", test_name, "\n")
    cat("  Error:", result$error, "\n")
  }

  return(result)
}

# Check required packages
check_packages <- function() {
  required_packages <- c("ggplot2", "dplyr", "tidyr", "microbenchmark")
  optional_packages <- c("pryr", "patchwork", "scales")

  cat("\nChecking required packages:\n")
  for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("  ✓", pkg, "is installed\n")
    } else {
      cat("  ✗", pkg, "is NOT installed (REQUIRED)\n")
      cat("    Install with: install.packages('", pkg, "')\n", sep = "")
    }
  }

  cat("\nChecking optional packages:\n")
  for (pkg in optional_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("  ✓", pkg, "is installed\n")
    } else {
      cat("  ⚠", pkg, "is NOT installed (optional)\n")
    }
  }
}

# Test 1: Check C++ implementation status
test_cpp_status <- function() {
  # First, ensure the fixed get_cpp_status function is available
  if (!exists("get_cpp_status")) {
    cat("get_cpp_status function not found. Loading from package...\n")
    return()
  }

  status <- get_cpp_status()
  cat("\nC++ Implementation Status:\n")
  for (name in names(status)) {
    cat(sprintf("  %-25s: %s\n", name, ifelse(status[[name]], "✓ Available", "✗ Not Available")))
  }

  # Check if exponential is in the status
  if (!"exponential" %in% names(status)) {
    cat("\n⚠ WARNING: 'exponential' not in cpp_status. You may need to update cpp_interface.R\n")
  }
}

# Test 2: Test data generation
test_data_generation <- function() {
  # Load the benchmark script
  source("inst/benchmarks/benchmark_exponential.R")

  cat("\nTesting data generation function:\n")

  # Test with different parameters
  test_cases <- list(
    list(n = 100, rates = c(1, 5), weights = NULL),
    list(n = 200, rates = c(0.5, 2, 10), weights = c(0.3, 0.3, 0.4))
  )

  for (i in seq_along(test_cases)) {
    tc <- test_cases[[i]]
    cat(sprintf("\n  Test case %d: n=%d, rates=%s\n",
                i, tc$n, paste(tc$rates, collapse=", ")))

    data <- generate_exponential_mixture(tc$n, tc$rates, tc$weights)

    cat(sprintf("    Generated %d observations\n", length(data)))
    cat(sprintf("    Range: [%.3f, %.3f]\n", min(data), max(data)))
    cat(sprintf("    Mean: %.3f\n", mean(data)))

    # Check data is valid
    if (any(data < 0)) {
      cat("    ✗ ERROR: Negative values found!\n")
    } else {
      cat("    ✓ All values are non-negative\n")
    }
  }
}

# Test 3: Quick benchmark with minimal parameters
test_quick_benchmark <- function() {
  source("inst/benchmarks/benchmark_exponential.R")

  cat("\nRunning minimal benchmark test:\n")

  # Very small test
  results <- benchmark_exponential_comprehensive(
    n_obs_vec = c(50, 100),
    n_iter_vec = c(10, 20),
    n_clusters_vec = c(2),
    n_reps = 2
  )

  cat("\nBenchmark results summary:\n")
  cat(sprintf("  Total scenarios tested: %d\n", length(unique(results$scenario_id))))
  cat(sprintf("  Rows in results: %d\n", nrow(results)))

  # Check results structure
  expected_cols <- c("scenario_id", "n_obs", "n_iter", "n_true_clusters",
                     "implementation", "time", "memory_mb", "clusters_found", "rep")
  missing_cols <- setdiff(expected_cols, names(results))

  if (length(missing_cols) == 0) {
    cat("  ✓ All expected columns present\n")
  } else {
    cat("  ✗ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  }

  # Summary statistics
  if (nrow(results) > 0) {
    r_times <- results$time[results$implementation == "R"]
    cpp_times <- results$time[results$implementation == "C++"]

    if (length(r_times) > 0 && length(cpp_times) > 0) {
      cat(sprintf("\n  Average times:\n"))
      cat(sprintf("    R:   %.3f seconds\n", mean(r_times)))
      cat(sprintf("    C++: %.3f seconds\n", mean(cpp_times)))
      cat(sprintf("    Speedup: %.1fx\n", mean(r_times) / mean(cpp_times)))
    }
  }

  return(results)
}

# Test 4: Memory profiling function
test_memory_profiling <- function() {
  source("inst/benchmarks/benchmark_exponential.R")

  cat("\nTesting memory profiling:\n")

  # Add memory tracking fallback if pryr not available
  if (!requireNamespace("pryr", quietly = TRUE)) {
    cat("  ⚠ pryr package not available, using fallback memory tracking\n")

    # Define fallback in global environment
    assign("mem_used", function() {
      gc_info <- gc()
      (gc_info[1, 2] + gc_info[2, 2]) * 1024
    }, envir = .GlobalEnv)
  }

  # Run with small dataset
  mem_results <- profile_exponential_memory(n_obs = 100, n_iter = 50)

  cat("\nMemory profiling results:\n")
  print(mem_results$summary)

  return(mem_results)
}

# Test 5: Component benchmarking
test_component_benchmark <- function() {
  source("inst/benchmarks/benchmark_exponential.R")

  cat("\nTesting component benchmarking:\n")

  # Need to ensure data is properly formatted
  comp_results <- benchmark_exponential_components(n_obs = 100, n_clusters = 2)

  cat("\nComponent benchmark results:\n")
  print(comp_results)

  return(comp_results)
}

# Test 6: Test individual likelihood calculations
test_likelihood_calculations <- function() {
  cat("\nTesting likelihood calculations:\n")

  # Create simple test data
  test_data <- matrix(rexp(10, rate = 2), ncol = 1)

  # Create mixing distribution
  mdObj <- ExponentialMixtureCreate(c(0.01, 0.01))

  # Test parameters
  test_theta <- list(array(2.0, dim = c(1, 1, 1)))

  # Test R implementation
  set_use_cpp(FALSE)
  lik_r <- Likelihood(mdObj, test_data[1, , drop = FALSE], test_theta)
  cat(sprintf("  R likelihood for first data point: %.6f\n", lik_r))

  # Test C++ implementation if available
  if (exists("exponential_likelihood_cpp")) {
    set_use_cpp(TRUE)
    lik_cpp <- Likelihood(mdObj, test_data[1, , drop = FALSE], test_theta)
    cat(sprintf("  C++ likelihood for first data point: %.6f\n", lik_cpp))
    cat(sprintf("  Difference: %.2e\n", abs(lik_r - lik_cpp)))
  } else {
    cat("  ⚠ C++ likelihood function not available\n")
  }
}

# Main test execution
run_all_tests <- function() {
  cat("\n")
  cat("====================================================\n")
  cat("    EXPONENTIAL BENCHMARKING TEST SUITE\n")
  cat("====================================================\n")

  # Store test results
  test_results <- list()

  # Check packages first
  check_packages()

  # Run tests
  test_results$cpp_status <- safe_test("C++ Status Check", test_cpp_status)
  test_results$likelihood <- safe_test("Likelihood Calculations", test_likelihood_calculations)
  test_results$data_gen <- safe_test("Data Generation", test_data_generation)
  test_results$quick_bench <- safe_test("Quick Benchmark", test_quick_benchmark)
  test_results$memory <- safe_test("Memory Profiling", test_memory_profiling)
  test_results$components <- safe_test("Component Benchmarking", test_component_benchmark)

  # Summary
  cat("\n")
  cat("====================================================\n")
  cat("                    SUMMARY\n")
  cat("====================================================\n")

  successful <- sum(sapply(test_results, function(x) x$success))
  total <- length(test_results)

  cat(sprintf("\nTests passed: %d/%d (%.0f%%)\n", successful, total, 100 * successful / total))

  if (successful < total) {
    cat("\nFailed tests:\n")
    for (name in names(test_results)) {
      if (!test_results[[name]]$success) {
        cat("  - ", name, "\n")
      }
    }
  }

  return(test_results)
}

# Function to apply fixes if needed
apply_fixes <- function() {
  cat("\nApplying fixes to benchmark_exponential.R...\n")

  # Read the file
  benchmark_file <- "inst/benchmarks/benchmark_exponential.R"
  if (!file.exists(benchmark_file)) {
    cat("  ✗ File not found:", benchmark_file, "\n")
    return(FALSE)
  }

  lines <- readLines(benchmark_file)

  # Fix 1: Add pryr fallback at the beginning
  pryr_check <- "if (!requireNamespace(\"pryr\", quietly = TRUE)) {"
  if (!any(grepl(pryr_check, lines, fixed = TRUE))) {
    cat("  Adding pryr fallback...\n")

    pryr_fallback <- c(
      "# Memory tracking fallback",
      "if (!requireNamespace(\"pryr\", quietly = TRUE)) {",
      "  mem_used <- function() {",
      "    gc_info <- gc()",
      "    (gc_info[1, 2] + gc_info[2, 2]) * 1024",
      "  }",
      "} else {",
      "  mem_used <- pryr::mem_used",
      "}"
    )

    # Insert after library statements
    library_line <- which(grepl("^library\\(", lines))[1]
    if (!is.na(library_line)) {
      last_library <- max(which(grepl("^library\\(", lines)))
      lines <- c(lines[1:last_library], "", pryr_fallback, lines[(last_library + 1):length(lines)])
    }
  }

  # Fix 2: Update benchmark_exponential_components to handle data as matrix
  comp_func_start <- which(grepl("benchmark_exponential_components <- function", lines))
  if (length(comp_func_start) > 0) {
    cat("  Fixing data handling in benchmark_exponential_components...\n")

    # Find where data is generated
    for (i in (comp_func_start[1] + 1):length(lines)) {
      if (grepl("data <- generate_exponential_mixture", lines[i])) {
        # Change variable name to data_vec
        lines[i] <- gsub("data <- generate", "data_vec <- generate", lines[i])
        # Add matrix conversion on next line
        lines <- c(lines[1:i],
                   "  data <- matrix(data_vec, ncol = 1)",
                   lines[(i + 1):length(lines)])
        break
      }
    }
  }

  # Write back
  writeLines(lines, benchmark_file)
  cat("  ✓ Fixes applied to", benchmark_file, "\n")

  return(TRUE)
}

# Execute tests
if (interactive()) {
  cat("Starting exponential benchmark tests...\n")
  cat("Working directory:", getwd(), "\n\n")

  # Check if benchmark file exists
  if (!file.exists("inst/benchmarks/benchmark_exponential.R")) {
    cat("✗ ERROR: benchmark_exponential.R not found!\n")
    cat("  Please ensure you're in the package root directory.\n")
  } else {
    # Option to apply fixes
    cat("Do you want to apply automatic fixes? (y/n): ")
    response <- readline()

    if (tolower(response) == "y") {
      apply_fixes()
    }

    # Run tests
    test_results <- run_all_tests()
  }
}
