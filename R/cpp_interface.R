#' @title C++ Backend Interface Functions
#' @description Functions for interfacing with C++ implementations
#' @name cpp_interface
NULL

#' Set whether to use C++ implementations
#' @param use_cpp Logical indicating whether to use C++ implementations
#' @export
set_use_cpp <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp = use_cpp)
  invisible(use_cpp)
}

#' Check if using C++ implementations
#' @return Logical indicating whether C++ implementations are being used
#' @export
using_cpp <- function() {
  getOption("dirichletprocess.use_cpp", FALSE)
}

#' Get C++ implementation status
#' @return List showing which C++ implementations are available
#' @export
get_cpp_status <- function() {
  status <- list(
    mcmc_runner = exists("_dirichletprocess_run_mcmc_cpp"),
    gaussian_likelihood = exists("_dirichletprocess_run_mcmc_cpp"),
    exponential_likelihood = exists("_dirichletprocess_run_mcmc_cpp"),
    beta_likelihood = exists("_dirichletprocess_run_mcmc_cpp"),
    hierarchical_beta = exists("_dirichletprocess_hierarchical_beta_fit_cpp"),
    benchmark_components = exists("_dirichletprocess_benchmark_components_cpp"),
    memory_tracking = exists("_dirichletprocess_get_memory_tracking")
  )
  return(status)
}

#' Check if C++ can be used for a given DP object
#' @param dp_obj Dirichlet process object
#' @return Logical indicating whether C++ implementation is available
#' @keywords internal
can_use_cpp <- function(dp_obj) {
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    return(FALSE)
  }

  supported_types <- c("normal_inverse_gamma", "normal", "mvnormal")
  inherits(dp_obj$mixingDistribution, supported_types)
}

#' Enable C++ implementations for specific samplers
#' @export
enable_cpp_samplers <- function() {
  invisible(exists("_dirichletprocess_run_mcmc_cpp", mode = "function"))
}

#' Check if using C++ samplers
#' @export
using_cpp_samplers <- function() {
  using_cpp() && exists("_dirichletprocess_run_mcmc_cpp", mode = "function")
}

#' Enable C++ implementations for hierarchical models
#' @export
enable_cpp_hierarchical_samplers <- function() {
  invisible(exists("_dirichletprocess_hierarchical_beta_fit_cpp", mode = "function"))
}

#' Check if using C++ hierarchical samplers
#' @export
using_cpp_hierarchical_samplers <- function() {
  using_cpp() && exists("_dirichletprocess_hierarchical_beta_fit_cpp", mode = "function")
}

#' Enable C++ implementations for Markov models
#' @export
enable_cpp_markov_samplers <- function() {
  invisible(exists("_dirichletprocess_markov_dp_fit_cpp", mode = "function"))
}

#' Check if using C++ Markov samplers
#' @export
using_cpp_markov_samplers <- function() {
  using_cpp() && exists("_dirichletprocess_markov_dp_fit_cpp", mode = "function")
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
  result$likelihoodChain <- as.numeric(result$likelihood_chain)

  return(result)
}

#' Create mixing distribution parameters for C++
#' @keywords internal
# Add to prepare_mixing_dist_params function
prepare_mixing_dist_params <- function(dp_obj) {
  md <- dp_obj$mixingDistribution

  if (inherits(md, "mvnormal")) {
    # Extract parameters from the mixing distribution
    if (!is.null(md$priorParameters)) {
      pp <- md$priorParameters
      list(
        type = "mvnormal",
        mu0 = as.numeric(pp$mu0),
        kappa0 = as.numeric(pp$kappa0),
        Lambda = as.matrix(pp$Lambda),
        nu = as.numeric(pp$nu)
      )
    } else {
      stop("MVNormal mixing distribution missing prior parameters")
    }
  } else if (inherits(md, "beta")) {
    # Beta distribution parameters
    list(
      type = "beta",
      alpha0 = md$priorParameters[1],
      beta0 = md$priorParameters[2],
      maxT = ifelse(is.null(md$maxT), 1, md$maxT),
      mhStepSize = md$mhStepSize,
      hyperPriorParameters = md$hyperPriorParameters
    )
  } else if (inherits(md, "normal_inverse_gamma") || inherits(md, "normal")) {
    # Existing Gaussian implementation
    if (!is.null(md$priors)) {
      list(
        type = "gaussian",
        mu0 = md$priors$mu_0,
        kappa0 = md$priors$kappa_0,
        alpha0 = md$priors$alpha_0,
        beta0 = md$priors$beta_0
      )
    } else {
      list(
        type = "gaussian",
        mu0 = md$priorParameters[1],
        kappa0 = md$priorParameters[2],
        alpha0 = md$priorParameters[3],
        beta0 = md$priorParameters[4]
      )
    }
  } else if (inherits(md, "exponential")) {
    list(
      type = "exponential",
      alpha0 = md$priorParameters[1],
      beta0 = md$priorParameters[2]
    )
  } else {
    stop("Mixing distribution not yet implemented in C++: ", class(md))
  }
}

#' Create MCMC parameters for C++
#' @keywords internal
prepare_mcmc_params <- function(dp_obj, its, updatePrior, n_burn = 0, thin = 1) {
  mcmc_params <- list(
    n_iter = as.integer(its),
    n_burn = as.integer(n_burn),
    thin = as.integer(thin),
    update_concentration = as.logical(updatePrior),
    alpha = as.numeric(dp_obj$alpha),
    m_auxiliary = 3L  # Default for Algorithm 8
  )

  # Add alpha prior parameters
  if (!is.null(dp_obj$alphaPriorParameters)) {
    mcmc_params$alpha_prior_shape <- dp_obj$alphaPriorParameters[1]
    mcmc_params$alpha_prior_rate <- dp_obj$alphaPriorParameters[2]
  } else {
    # Default Gamma(1,1) prior
    mcmc_params$alpha_prior_shape <- 1.0
    mcmc_params$alpha_prior_rate <- 1.0
  }

  return(mcmc_params)
}
