# tests/testthat/helper-debug.R
# Debug utilities for MCMC testing

#' Comprehensive system check for C++ availability
check_cpp_system <- function() {
  cat("\n=== SYSTEM CHECK ===\n")

  # Basic checks
  cat("R version:", R.version.string, "\n")
  cat("Platform:", .Platform$OS.type, "\n")

  # Package checks
  if (requireNamespace("dirichletprocess", quietly = TRUE)) {
    cat("dirichletprocess package loaded: YES\n")
  } else {
    cat("dirichletprocess package loaded: NO\n")
    return(FALSE)
  }

  # C++ function checks
  cpp_functions <- c(
    "_dirichletprocess_run_mcmc_cpp",
    "_dirichletprocess_normal_prior_draw_cpp",
    "_dirichletprocess_normal_posterior_draw_cpp"
  )

  for (func in cpp_functions) {
    cat("Function", func, "exists:", exists(func), "\n")
  }

  # Backend status
  cpp_status <- get_cpp_status()
  cat("get_cpp_status():", paste(names(cpp_status), cpp_status, sep = "=", collapse = ", "), "\n")

  # Test simple C++ call
  tryCatch({
    set_use_cpp(TRUE)
    simple_dp <- DirichletProcessGaussian(c(1, 2, 3))
    cat("Simple C++ DP creation: SUCCESS\n")
    cat("can_use_cpp():", can_use_cpp(simple_dp), "\n")
  }, error = function(e) {
    cat("Simple C++ DP creation: FAILED -", e$message, "\n")
  })

  return(TRUE)
}

#' Safe value printing for debug
safe_print_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  } else if (is.list(x)) {
    return(paste0("LIST(", length(x), " elements)"))
  } else if (is.numeric(x) && length(x) == 1) {
    return(as.character(x))
  } else if (is.numeric(x) && length(x) > 1) {
    return(paste0("VECTOR(", length(x), " elements)"))
  } else {
    return(paste0("OTHER(", class(x)[1], ")"))
  }
}

#' Deep inspection of DP objects - FIXED VERSION
inspect_dp_object <- function(dp_obj, name = "DP") {
  cat("\n=== INSPECTING", name, "===\n")
  cat("Class:", paste(class(dp_obj), collapse = " "), "\n")
  cat("Names:", paste(names(dp_obj), collapse = ", "), "\n")
  cat("Number of clusters:", dp_obj$numberClusters %||% dp_obj$n_clusters %||% "NULL", "\n")

  if ("alpha" %in% names(dp_obj)) {
    cat("Alpha value:", safe_print_value(dp_obj$alpha), "class:", paste(class(dp_obj$alpha), collapse = " "), "\n")
  }

  if ("alphaChain" %in% names(dp_obj)) {
    cat("AlphaChain length:", length(dp_obj$alphaChain), "\n")
    if (length(dp_obj$alphaChain) > 0) {
      cat("AlphaChain range:", paste(range(dp_obj$alphaChain), collapse = " to "), "\n")
    }
  } else {
    cat("AlphaChain: NOT PRESENT\n")
  }

  if ("clusterLabels" %in% names(dp_obj)) {
    cat("ClusterLabels length:", length(dp_obj$clusterLabels), "\n")
    cat("Unique labels:", length(unique(dp_obj$clusterLabels)), "\n")
  }

  if ("cluster_labels" %in% names(dp_obj)) {
    cat("cluster_labels (C++ style) length:", length(dp_obj$cluster_labels), "\n")
  }

  if ("likelihoodChain" %in% names(dp_obj)) {
    cat("LikelihoodChain length:", length(dp_obj$likelihoodChain), "\n")
  } else {
    cat("LikelihoodChain: NOT PRESENT\n")
  }

  # Check for C++ specific fields
  cpp_fields <- c("theta", "n_clusters", "cluster_labels")
  present_cpp <- cpp_fields[cpp_fields %in% names(dp_obj)]
  if (length(present_cpp) > 0) {
    cat("C++ specific fields present:", paste(present_cpp, collapse = ", "), "\n")
  }
}

#' Test data generator with validation
create_test_data <- function(n = 50, k = 2, sep = 3, seed = 123) {
  set.seed(seed)

  # Create clearly separated clusters
  cluster_means <- seq(-sep, sep, length.out = k)
  cluster_sizes <- rmultinom(1, n, rep(1/k, k))[, 1]

  data <- c()
  labels <- c()

  for (i in 1:k) {
    if (cluster_sizes[i] > 0) {
      cluster_data <- rnorm(cluster_sizes[i], mean = cluster_means[i], sd = 0.5)
      data <- c(data, cluster_data)
      labels <- c(labels, rep(i, cluster_sizes[i]))
    }
  }

  cat("Generated test data: n =", length(data), "true k =", k,
      "separation =", sep, "\n")
  cat("True cluster sizes:", paste(table(labels), collapse = " "), "\n")

  return(list(data = data, true_labels = labels, true_k = k))
}

#' Safe implementation comparison
safe_compare_implementations <- function(data, n_iter = 50, seed = 123) {
  results <- list()

  # Always try R first
  tryCatch({
    set_use_cpp(FALSE)
    set.seed(seed)
    dp_r <- DirichletProcessGaussian(data)
    dp_r <- Fit(dp_r, n_iter, progressBar = FALSE)
    results$r <- dp_r
    results$r_success <- TRUE
    cat("R implementation: SUCCESS\n")
  }, error = function(e) {
    cat("R implementation: FAILED -", e$message, "\n")
    results$r_success <- FALSE
  })

  # Try C++ if available
  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    tryCatch({
      set_use_cpp(TRUE)
      set.seed(seed)
      dp_cpp <- DirichletProcessGaussian(data)
      dp_cpp <- Fit(dp_cpp, n_iter, progressBar = FALSE)
      results$cpp <- dp_cpp
      results$cpp_success <- TRUE
      cat("C++ implementation: SUCCESS\n")
    }, error = function(e) {
      cat("C++ implementation: FAILED -", e$message, "\n")
      results$cpp_success <- FALSE
    })
  } else {
    results$cpp_success <- FALSE
    cat("C++ implementation: NOT AVAILABLE\n")
  }

  results$both_available <- results$r_success && results$cpp_success
  return(results)
}
