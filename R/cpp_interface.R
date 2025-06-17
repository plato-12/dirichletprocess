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

#' Run MCMC using C++ implementation with proper chain handling
#' @keywords internal
run_mcmc_cpp <- function(data, mixing_dist_params, mcmc_params) {
  # Ensure mcmc_params has all required fields
  if (!"m_auxiliary" %in% names(mcmc_params)) {
    mcmc_params$m_auxiliary <- 3  # Default value matching R implementation
  }

  if (!"alpha" %in% names(mcmc_params)) {
    mcmc_params$alpha <- 1.0  # Default starting value
  }

  if (!"update_concentration" %in% names(mcmc_params)) {
    mcmc_params$update_concentration <- TRUE
  }

  # Call C++ implementation
  result <- .Call("_dirichletprocess_run_mcmc_cpp",
                  data, mixing_dist_params, mcmc_params,
                  PACKAGE = "dirichletprocess")

  # Convert results to match R format
  result$labelsChain <- lapply(1:nrow(result$labels_chain), function(i) {
    result$labels_chain[i,]
  })

  result$alphaChain <- as.numeric(result$alpha_chain)

  return(result)
}

#' Create mixing distribution parameters for C++
#' @keywords internal
prepare_mixing_dist_params <- function(dp_obj) {
  md <- dp_obj$mixingDistribution

  # Handle normal/gaussian distributions
  if (inherits(md, c("normal", "normal_inverse_gamma", "gaussian", "conjugate"))) {
    # Check for different parameter structures
    if (!is.null(md$priors)) {
      # New structure with named priors
      return(list(
        type = "gaussian",
        mu0 = as.numeric(md$priors$mu_0),
        kappa0 = as.numeric(md$priors$kappa_0),
        alpha0 = as.numeric(md$priors$alpha_0),
        beta0 = as.numeric(md$priors$beta_0)
      ))
    } else if (!is.null(md$priorParameters)) {
      # Old structure with priorParameters vector
      params <- md$priorParameters
      if (is.list(params)) {
        # If it's a list, extract values
        return(list(
          type = "gaussian",
          mu0 = as.numeric(params[[1]]),
          kappa0 = as.numeric(params[[2]]),
          alpha0 = as.numeric(params[[3]]),
          beta0 = as.numeric(params[[4]])
        ))
      } else if (is.numeric(params) && length(params) == 4) {
        # If it's a numeric vector
        return(list(
          type = "gaussian",
          mu0 = as.numeric(params[1]),
          kappa0 = as.numeric(params[2]),
          alpha0 = as.numeric(params[3]),
          beta0 = as.numeric(params[4])
        ))
      }
    }

    # Default parameters if none found
    warning("Using default prior parameters for normal distribution")
    return(list(
      type = "gaussian",
      mu0 = 0.0,
      kappa0 = 1.0,
      alpha0 = 1.0,
      beta0 = 1.0
    ))
  }

  stop("Mixing distribution not yet implemented in C++: ", class(md))
}

#' Check if C++ implementation is available for this model
#' @keywords internal
can_use_cpp <- function(dp_obj) {
  # Check if C++ function exists in multiple ways
  cpp_available <- exists("run_mcmc_cpp", where = asNamespace("dirichletprocess"),
                          mode = "function") ||
    exists("_dirichletprocess_run_mcmc_cpp",
           where = asNamespace("dirichletprocess")) ||
    exists("_dirichletprocess_run_mcmc_cpp")

  if (!cpp_available) {
    return(FALSE)
  }

  # Check if the mixing distribution is supported
  supported_types <- c("normal", "normal_inverse_gamma", "gaussian", "conjugate")

  if (inherits(dp_obj$mixingDistribution, supported_types)) {
    return(TRUE)
  }

  # Additional check for normal distribution with different class structures
  if (!is.null(dp_obj$mixingDistribution$priorParameters) &&
      length(dp_obj$mixingDistribution$priorParameters) == 4) {
    return(TRUE)
  }

  return(FALSE)
}

prepare_mcmc_params <- function(n_iter, n_burn = 0, thin = 1,
                                update_concentration = TRUE,
                                alpha = NULL, ...) {
  # If alpha not specified, use a reasonable default based on data size
  if (is.null(alpha)) {
    # Use a default that encourages exploration
    alpha <- 1.0
  }

  list(
    n_iter = as.integer(n_iter),
    n_burn = as.integer(n_burn),
    thin = as.integer(thin),
    update_concentration = as.logical(update_concentration),
    alpha = as.numeric(alpha)
  )
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
