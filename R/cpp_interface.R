#' Set whether to use C++ implementations
#'
#' @param use_cpp Logical indicating whether to use C++ implementations
#' @export
set_use_cpp <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp = use_cpp)

  if (use_cpp) {
    message("C++ implementations enabled")
  } else {
    message("Using R implementations")
  }

  invisible(use_cpp)
}

#' Check if C++ implementations are enabled
#'
#' @return Logical indicating if C++ implementations are enabled
#' @export
using_cpp <- function() {
  getOption("dirichletprocess.use_cpp", FALSE)
}

#' Get C++ Implementation Status
#'
#' Check which C++ implementations are available
#'
#' @return A named logical vector indicating which C++ implementations are available
#' @export
get_cpp_status <- function() {
  # Check for existence of key C++ functions
  status <- list(
    # Core MCMC runner - check for the actual Rcpp export
    mcmc_runner = tryCatch({
      exists("_dirichletprocess_run_mcmc_cpp") ||
        exists(".Call") && exists("run_mcmc_cpp", envir = getNamespace("dirichletprocess"))
    }, error = function(e) FALSE),

    # Distribution-specific functions
    gaussian_likelihood = tryCatch({
      exists("_dirichletprocess_run_mcmc_cpp")
    }, error = function(e) FALSE),

    # Legacy functions for backward compatibility
    likelihood = exists("likelihood_cpp", mode = "function"),
    normal_likelihood = exists("normal_likelihood_cpp", mode = "function"),
    exponential_likelihood = exists("exponential_likelihood_cpp", mode = "function"),
    mvnormal_likelihood = exists("mvnormal_likelihood_cpp", mode = "function"),

    # Hierarchical models
    hierarchical_beta = exists("hierarchical_beta_fit_cpp", mode = "function"),
    hierarchical_mvnormal = exists("hierarchical_mvnormal2_fit_cpp", mode = "function"),

    # Markov DP
    markov_dp = exists("markov_dp_fit_cpp", mode = "function"),

    # Benchmarking and profiling
    benchmark_components = exists("benchmark_cpp_components", mode = "function"),
    memory_tracking = exists("get_memory_tracking", mode = "function")
  )

  return(status)
}

#' Check if C++ Implementation is Available
#'
#' @param func_name Name of the function to check
#' @return Logical indicating if the C++ implementation is available
#' @keywords internal
has_cpp_implementation <- function(func_name) {
  exists(func_name, mode = "function") ||
    exists(paste0("_dirichletprocess_", func_name), mode = "function")
}

#' Initialize C++ availability check
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Set default options
  options(dirichletprocess.use_cpp = FALSE)
  options(dirichletprocess.verbose = FALSE)
}

#' Package startup message
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Check if C++ implementations are available
  cpp_status <- tryCatch({
    get_cpp_status()
  }, error = function(e) {
    list(mcmc_runner = FALSE)
  })

  available_count <- sum(unlist(cpp_status))

  if (available_count > 0) {
    if (getOption("dirichletprocess.verbose", FALSE)) {
      packageStartupMessage(
        sprintf("dirichletprocess: %d C++ implementations available", available_count)
      )
      packageStartupMessage("Use set_use_cpp(TRUE) to enable C++ implementations")
    }
  } else {
    if (getOption("dirichletprocess.verbose", FALSE)) {
      packageStartupMessage("dirichletprocess: Using R implementations (C++ not available)")
    }
  }
}

#' Run MCMC using C++ implementation
#' @keywords internal
run_mcmc_cpp <- function(data, mixing_dist_params, mcmc_params) {
  # FIXED: Check for the actual Rcpp exported function name
  if (exists("_dirichletprocess_run_mcmc_cpp")) {
    return(.Call("_dirichletprocess_run_mcmc_cpp",
                 data = as.matrix(data),
                 mixing_dist_params = mixing_dist_params,
                 mcmc_params = mcmc_params,
                 PACKAGE = "dirichletprocess"))
  } else {
    stop("C++ MCMC implementation not available. Make sure the package is compiled with C++ support.")
  }
}

#' Create mixing distribution parameters for C++
#' @keywords internal
prepare_mixing_dist_params <- function(dp_obj) {
  md <- dp_obj$mixingDistribution

  if (inherits(md, "normal_inverse_gamma") || inherits(md, "normal")) {
    # Handle both old and new parameter structures
    if (!is.null(md$priors)) {
      # New structure
      list(
        type = "gaussian",
        mu0 = md$priors$mu_0,
        kappa0 = md$priors$kappa_0,
        alpha0 = md$priors$alpha_0,
        beta0 = md$priors$beta_0
      )
    } else {
      # Old structure with priorParameters
      list(
        type = "gaussian",
        mu0 = md$priorParameters[1],
        kappa0 = md$priorParameters[2],
        alpha0 = md$priorParameters[3],
        beta0 = md$priorParameters[4]
      )
    }
  } else {
    stop("Mixing distribution not yet implemented in C++: ", class(md))
  }
}

#' Check if C++ implementation is available for this model
#' @keywords internal
can_use_cpp <- function(dp_obj) {
  # For now, only return TRUE if we actually have the C++ implementation
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    return(FALSE)
  }

  supported_types <- c("normal_inverse_gamma", "normal")
  inherits(dp_obj$mixingDistribution, supported_types)
}

#' Enable/disable cpp hierarchical samplers
#' @keywords internal
enable_cpp_hierarchical_samplers <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp_hierarchical = use_cpp)
  invisible(use_cpp)
}

#' Check if cpp hierarchical samplers are enabled
#' @keywords internal
using_cpp_hierarchical_samplers <- function() {
  getOption("dirichletprocess.use_cpp_hierarchical", FALSE)
}

#' Enable/disable cpp markov samplers
#' @keywords internal
enable_cpp_markov_samplers <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp_markov = use_cpp)
  invisible(use_cpp)
}

#' Check if cpp markov samplers are enabled
#' @keywords internal
using_cpp_markov_samplers <- function() {
  getOption("dirichletprocess.use_cpp_markov", FALSE)
}

#' Enable/disable cpp samplers (general)
#' @keywords internal
enable_cpp_samplers <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp_samplers = use_cpp)
  invisible(use_cpp)
}

#' Check if cpp samplers are enabled (general)
#' @keywords internal
using_cpp_samplers <- function() {
  getOption("dirichletprocess.use_cpp_samplers", FALSE)
}
