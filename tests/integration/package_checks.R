# tests/integration/package_checks.R

run_package_checks <- function() {
  cat("\n=== RUNNING PACKAGE CHECKS ===\n")

  check_results <- list()

  # 1. Run devtools::test()
  cat("\n1. Running unit tests...\n")
  test_result <- devtools::test()
  check_results$tests <- test_result

  # 2. Run devtools::check()
  cat("\n2. Running R CMD check...\n")
  check_result <- devtools::check()
  check_results$check <- check_result

  # 3. Check for compilation warnings
  cat("\n3. Checking C++ compilation...\n")
  cpp_check <- check_cpp_compilation()
  check_results$cpp <- cpp_check

  # 4. Documentation check
  cat("\n4. Checking documentation...\n")
  doc_check <- devtools::document()
  check_results$documentation <- doc_check

  # 5. Example check
  cat("\n5. Running examples...\n")
  example_check <- devtools::run_examples()
  check_results$examples <- example_check

  # Create summary report
  create_check_summary(check_results)

  return(check_results)
}

check_cpp_compilation <- function() {
  # Clean and rebuild
  devtools::clean_dll()

  # Capture compilation output
  compile_log <- capture.output({
    devtools::compile_dll()
  })

  # Check for warnings
  warnings <- grep("warning:", compile_log, value = TRUE)
  errors <- grep("error:", compile_log, value = TRUE)

  list(
    success = length(errors) == 0,
    warnings = warnings,
    errors = errors,
    full_log = compile_log
  )
}

create_check_summary <- function(check_results) {
  cat("\n\n=== CHECK SUMMARY ===\n")

  # Test summary
  test_summary <- check_results$tests
  cat("\nUnit Tests:",
      test_summary$failed, "failures,",
      test_summary$warnings, "warnings,",
      test_summary$skipped, "skipped,",
      test_summary$passed, "passed\n")

  # R CMD check summary
  check_summary <- check_results$check
  cat("\nR CMD check:",
      length(check_summary$errors), "errors,",
      length(check_summary$warnings), "warnings,",
      length(check_summary$notes), "notes\n")

  # C++ compilation
  if (check_results$cpp$success) {
    cat("\nC++ Compilation: SUCCESS")
    if (length(check_results$cpp$warnings) > 0) {
      cat(" (", length(check_results$cpp$warnings), " warnings)")
    }
    cat("\n")
  } else {
    cat("\nC++ Compilation: FAILED\n")
  }

  # Overall status
  overall_success <-
    test_summary$failed == 0 &&
    length(check_summary$errors) == 0 &&
    check_results$cpp$success

  cat("\n", ifelse(overall_success, "✅ READY FOR RELEASE", "❌ ISSUES NEED FIXING"), "\n")
}

# Check specific functionality
check_cpp_availability <- function() {
  cat("\nChecking C++ availability for all distributions...\n")

  distributions <- list(
    normal = DirichletProcessGaussian,
    exponential = DirichletProcessExponential,
    beta = DirichletProcessBeta,
    weibull = DirichletProcessWeibull,
    mvnormal = DirichletProcessMvnormal,
    mvnormal2 = DirichletProcessMvnormal2
  )

  results <- list()

  for (dist_name in names(distributions)) {
    test_data <- generate_test_data(dist_name, n = 50)
    dp <- distributions[[dist_name]](test_data)

    # Check if C++ is available
    cpp_available <- can_use_cpp(dp)

    # Try to run with C++
    set_use_cpp(TRUE)
    cpp_works <- tryCatch({
      dp_test <- Fit(dp, its = 10)
      TRUE
    }, error = function(e) {
      FALSE
    })

    results[[dist_name]] <- list(
      available = cpp_available,
      works = cpp_works
    )

    cat("  ", dist_name, ": ",
        ifelse(cpp_available, "Available", "Not Available"),
        " | ",
        ifelse(cpp_works, "Works", "Fails"), "\n")
  }

  return(results)
}

# Check package dependencies
check_dependencies <- function() {
  cat("\nChecking package dependencies...\n")

  # Get package dependencies
  deps <- desc::desc_get_deps()

  # Check which are installed
  installed_pkgs <- installed.packages()[, "Package"]

  deps$installed <- deps$package %in% installed_pkgs

  # Check versions if installed
  for (i in 1:nrow(deps)) {
    if (deps$installed[i]) {
      deps$installed_version[i] <- as.character(packageVersion(deps$package[i]))
    }
  }

  print(deps)

  # Check for missing dependencies
  missing <- deps$package[!deps$installed]
  if (length(missing) > 0) {
    cat("\n⚠️  Missing dependencies:", paste(missing, collapse = ", "), "\n")
  } else {
    cat("\n✅ All dependencies installed\n")
  }

  return(deps)
}

# Validate NAMESPACE file
check_namespace <- function() {
  cat("\nChecking NAMESPACE file...\n")

  namespace_content <- readLines("NAMESPACE")

  # Check for C++ exports
  cpp_exports <- grep("^\\.Call|^useDynLib", namespace_content, value = TRUE)
  cat("C++ exports found:", length(cpp_exports), "\n")

  # Check for S3 methods
  s3_methods <- grep("^S3method", namespace_content, value = TRUE)
  cat("S3 methods registered:", length(s3_methods), "\n")

  # Check for imports
  imports <- grep("^import", namespace_content, value = TRUE)
  cat("Imports found:", length(imports), "\n")

  return(list(
    cpp_exports = cpp_exports,
    s3_methods = s3_methods,
    imports = imports
  ))
}

# Run all integration checks
run_all_integration_checks <- function() {
  results <- list()

  results$package_checks <- run_package_checks()
  results$cpp_availability <- check_cpp_availability()
  results$dependencies <- check_dependencies()
  results$namespace <- check_namespace()

  # Generate integration report
  generate_integration_report(results)

  return(results)
}

generate_integration_report <- function(results) {
  cat("\n\n=== INTEGRATION REPORT ===\n")
  cat("Generated:", format(Sys.time()), "\n\n")

  # Package check status
  cat("Package Checks:\n")
  cat("  R CMD check: ",
      ifelse(length(results$package_checks$check$errors) == 0, "✅ PASS", "❌ FAIL"), "\n")
  cat("  Unit tests: ",
      ifelse(results$package_checks$tests$failed == 0, "✅ PASS", "❌ FAIL"), "\n")
  cat("  C++ compilation: ",
      ifelse(results$package_checks$cpp$success, "✅ PASS", "❌ FAIL"), "\n")

  # C++ availability
  cat("\nC++ Implementation Status:\n")
  cpp_status <- results$cpp_availability
  for (dist in names(cpp_status)) {
    cat("  ", dist, ": ",
        ifelse(cpp_status[[dist]]$works, "✅", "❌"), "\n")
  }

  # Save detailed report
  saveRDS(results, "integration_results.rds")
  cat("\nDetailed results saved to: integration_results.rds\n")
}
