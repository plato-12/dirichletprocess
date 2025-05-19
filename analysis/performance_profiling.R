# Save as analysis/performance_profiling.R

# Attempt to ensure a clean environment for the package
if ("dirichletprocess" %in% loadedNamespaces()) {
  try(detach("package:dirichletprocess", unload = TRUE, character.only = TRUE, force = TRUE), silent = TRUE)
}
library(dirichletprocess)
library(profvis)
library(microbenchmark)

# Profile different mixture types with varying data sizes
profile_mixtures <- function() {
  # Ensure the results directory exists
  results_dir_prof <- "analysis/profiling_results" # Renamed for clarity
  if (!dir.exists(results_dir_prof)) {
    dir.create(results_dir_prof, recursive = TRUE)
  }

  profile_sizes <- c(100, 500, 1000)
  mcmc_iterations_profvis_conjugate <- 500
  mcmc_iterations_profvis_non_conjugate <- 200

  # Gaussian mixture (Conjugate)
  for (n in profile_sizes) {
    cat("Profiling Gaussian mixture with n =", n, ", Iterations =", mcmc_iterations_profvis_conjugate, "\n")
    set.seed(123)
    y_gauss <- rnorm(n) # Renamed for clarity

    p_gaussian <- NULL
    tryCatch({
      p_gaussian <- profvis({
        dp <- DirichletProcessGaussian(y_gauss)
        dp <- Fit(dp, mcmc_iterations_profvis_conjugate, progressBar = FALSE)
      })
      if (!is.null(p_gaussian)) {
        htmlwidgets::saveWidget(p_gaussian, file.path(results_dir_prof, paste0("gaussian_n", n, ".html")))
        cat("Saved Gaussian profiling for n =", n, "\n")
      } else {
        cat("Profvis object was NULL for Gaussian n =", n, ". Skipping save.\n")
      }
    }, error = function(e) {
      cat("Error profiling Gaussian for n =", n, ":", conditionMessage(e), "\n")
    })
  }

  # Beta mixture (Non-conjugate)
  for (n in profile_sizes[1:2]) {
    cat("Profiling Beta mixture with n =", n, ", Iterations =", mcmc_iterations_profvis_non_conjugate, "\n")
    set.seed(123)
    y_beta <- rbeta(n, 2, 5)

    p_beta <- NULL
    tryCatch({
      # Explicitly check if DirichletProcessBeta exists and its formals
      if(exists("DirichletProcessBeta", where="package:dirichletprocess")) {
        # cat("Formals for DirichletProcessBeta:", paste(names(formals(get("DirichletProcessBeta", envir = asNamespace("dirichletprocess")))), collapse=", "), "\n")
      }

      p_beta <- profvis({
        # Call using the package namespace explicitly IF simple call fails, though library() should suffice
        dp <- dirichletprocess::DirichletProcessBeta(y = y_beta, maxY = 1) # Pass y_beta as y
        dp <- Fit(dp, mcmc_iterations_profvis_non_conjugate, progressBar = FALSE)
      }, interval = 0.005) # Try smaller interval for profvis if functions are fast

      if (!is.null(p_beta)) {
        htmlwidgets::saveWidget(p_beta, file.path(results_dir_prof, paste0("beta_n", n, ".html")))
        cat("Saved Beta profiling for n =", n, "\n")
      } else {
        cat("Profvis object was NULL for Beta n =", n, ". Error likely occurred within profvis block.\n")
      }
    }, error = function(e) {
      cat("Error profiling Beta for n =", n, ":", conditionMessage(e), "\n")
      # Print traceback if error occurs in profvis
      if(exists("p_beta") && inherits(p_beta, "profvis")){
        print(p_beta) # This might show the error location if profvis captured it
      }
    })
  }

  # Multivariate normal (Conjugate)
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    # library(mvtnorm) # Already loaded by dirichletprocess if used internally or called here
    for (n in profile_sizes[1:2]) {
      cat("Profiling MVN mixture with n =", n, ", Iterations =", mcmc_iterations_profvis_conjugate, "\n")
      set.seed(123)
      y_mvn <- mvtnorm::rmvnorm(n, c(0, 0), diag(2))

      p_mvn <- NULL
      tryCatch({
        p_mvn <- profvis({
          dp <- DirichletProcessMvnormal(y_mvn)
          dp <- Fit(dp, mcmc_iterations_profvis_conjugate, progressBar = FALSE)
        })
        if(!is.null(p_mvn)){
          htmlwidgets::saveWidget(p_mvn, file.path(results_dir_prof, paste0("mvn_n", n, ".html")))
          cat("Saved MVN profiling for n =", n, "\n")
        } else {
          cat("Profvis object was NULL for MVN n =", n, ". Error likely occurred within profvis block.\n")
        }
      }, error = function(e) {
        cat("Error profiling MVN for n =", n, ":", conditionMessage(e), "\n")
      })
    }
  }
}

# Benchmark specific functions
benchmark_components <- function() {
  set.seed(123)
  results_dir_bench <- "analysis/profiling_results" # Renamed for clarity
  if (!dir.exists(results_dir_bench)) {
    dir.create(results_dir_bench, recursive = TRUE)
  }

  # --- Gaussian (Conjugate) Benchmarks ---
  cat("Benchmarking Gaussian components (n=200)\n")
  y_bench_gauss <- rnorm(200) # Renamed
  dp_gaussian_bench <- DirichletProcessGaussian(y_bench_gauss)
  dp_gaussian_bench <- Fit(dp_gaussian_bench, 5, progressBar=FALSE) # Increased init fit slightly


  bm1_gaussian <- NULL
  tryCatch({
    bm1_gaussian <- microbenchmark(
      ClusterComponentUpdate(dp_gaussian_bench),
      times = 20
    )
  }, error = function(e) {cat("Error benchmarking Gaussian ClusterComponentUpdate:", conditionMessage(e), "\n")})

  bm2_gaussian <- NULL
  tryCatch({
    bm2_gaussian <- microbenchmark(
      ClusterParameterUpdate(dp_gaussian_bench),
      times = 20
    )
  }, error = function(e) {cat("Error benchmarking Gaussian ClusterParameterUpdate:", conditionMessage(e), "\n")})

  bm3_gaussian <- NULL
  tryCatch({
    bm3_gaussian <- microbenchmark(
      UpdateAlpha(dp_gaussian_bench),
      times = 50
    )
  }, error = function(e) {cat("Error benchmarking Gaussian UpdateAlpha:", conditionMessage(e), "\n")})


  # --- Beta (Non-Conjugate) Benchmarks ---
  cat("\nBenchmarking Beta components (n=100 for potentially slower non-conjugate)\n")
  y_bench_beta <- rbeta(100, 2, 5) # Renamed
  dp_beta_bench <- NULL # Initialize

  tryCatch({
    # Crucial call that was failing
    dp_beta_bench <- dirichletprocess::DirichletProcessBeta(y = y_bench_beta, maxY = 1)
    dp_beta_bench <- Fit(dp_beta_bench, 5, progressBar=FALSE) # Increased init fit slightly
  }, error = function(e) {
    cat("CRITICAL ERROR initializing Beta DP for benchmarking:", conditionMessage(e), "\n")
    # dp_beta_bench will remain NULL, subsequent benchmarks for Beta will be skipped by !is.null checks
  })

  bm1_beta <- NULL
  if(!is.null(dp_beta_bench)) {
    tryCatch({
      bm1_beta <- microbenchmark(
        ClusterComponentUpdate(dp_beta_bench),
        times = 10
      )
    }, error = function(e) {cat("Error benchmarking Beta ClusterComponentUpdate:", conditionMessage(e), "\n")})
  }

  bm2_beta <- NULL
  if(!is.null(dp_beta_bench)) {
    tryCatch({
      bm2_beta <- microbenchmark(
        ClusterParameterUpdate(dp_beta_bench),
        times = 5
      )
    }, error = function(e) {cat("Error benchmarking Beta ClusterParameterUpdate:", conditionMessage(e), "\n")})
  }

  bm3_beta <- NULL
  if(!is.null(dp_beta_bench)) {
    tryCatch({
      bm3_beta <- microbenchmark(
        UpdateAlpha(dp_beta_bench),
        times = 50
      )
    }, error = function(e) {cat("Error benchmarking Beta UpdateAlpha:", conditionMessage(e), "\n")})
  }

  benchmark_results_list <- list(
    Gaussian_ClusterComponentUpdate = bm1_gaussian,
    Gaussian_ClusterParameterUpdate = bm2_gaussian,
    Gaussian_UpdateAlpha = bm3_gaussian
  )
  # Only add Beta results if dp_beta_bench was successfully created
  if(!is.null(dp_beta_bench)){
    benchmark_results_list$Beta_ClusterComponentUpdate = bm1_beta
    benchmark_results_list$Beta_ClusterParameterUpdate = bm2_beta
    benchmark_results_list$Beta_UpdateAlpha = bm3_beta
  }

  saveRDS(benchmark_results_list, file.path(results_dir_bench, "component_benchmarks.rds"))
  cat("\nBenchmark objects saved to analysis/profiling_results/component_benchmarks.rds\n")

  summary_file_path <- file.path(results_dir_bench, "benchmark_summary.txt")

  # Initialize or overwrite summary file
  file_conn <- file(summary_file_path, "w")
  writeLines("# Benchmark Results\n", file_conn)
  close(file_conn)

  append_benchmark_summary <- function(bm_obj, name, file_path_app) { # Renamed arg
    if (!is.null(bm_obj)) {
      sink(file_path_app, append = TRUE)
      cat("\n## ", name, "\n")
      print(summary(bm_obj)) # microbenchmark summary is a data.frame, print will format it
      cat("\n")
      sink()
    } else {
      sink(file_path_app, append = TRUE)
      cat("\n## ", name, "\n")
      cat("Benchmark object is NULL (likely due to an earlier error for this component).\n\n")
      sink()
    }
  }

  append_benchmark_summary(bm1_gaussian, "Gaussian ClusterComponentUpdate (n=200)", summary_file_path)
  append_benchmark_summary(bm2_gaussian, "Gaussian ClusterParameterUpdate (n=200)", summary_file_path)
  append_benchmark_summary(bm3_gaussian, "Gaussian UpdateAlpha (n=200)", summary_file_path)

  # Append Beta results only if benchmarks were run
  if(!is.null(dp_beta_bench)){
    append_benchmark_summary(bm1_beta, "Beta ClusterComponentUpdate (n=100, Non-conjugate)", summary_file_path)
    append_benchmark_summary(bm2_beta, "Beta ClusterParameterUpdate (n=100, Non-conjugate)", summary_file_path)
    append_benchmark_summary(bm3_beta, "Beta UpdateAlpha (n=100, Non-conjugate)", summary_file_path)
  } else {
    append_benchmark_summary(NULL, "Beta Benchmarks SKIPPED due to initialization error", summary_file_path)
  }


  cat("Benchmark summary written to analysis/profiling_results/benchmark_summary.txt\n")
}

# Run the profiling and benchmarking
profile_mixtures()
benchmark_components()

cat("Profiling and benchmarking complete.\nHTML profiles are in analysis/profiling_results/\nRDS benchmarks object and summary text file are also in analysis/profiling_results/\n")
