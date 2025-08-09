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
  has_cpp <- exists("_dirichletprocess_run_mcmc_cpp")

  status <- list(
    mcmc_runner = has_cpp,
    gaussian_likelihood = has_cpp,
    exponential_likelihood = exists("_dirichletprocess_run_mcmc_cpp"),
    beta_likelihood = has_cpp,
    mvnormal_likelihood = exists("conjugate_mvnormal_cluster_component_update_cpp"),
    weibull_likelihood = has_cpp,
    hierarchical_beta = exists("_dirichletprocess_hierarchical_beta_fit_cpp"),
    markov = exists("_dirichletprocess_markov_dp_fit_cpp"),
    available = has_cpp
  )

  attr(status, "message") <- if (has_cpp) {
    "C++ backend is available with all distributions"
  } else {
    "C++ backend is not available - using R implementation"
  }

  return(status)
}

#' Check if C++ can be used for a given DP object
#' @param dp_obj Dirichlet process object
#' @return Logical indicating whether C++ implementation is available
#' @export
can_use_cpp <- function(dp_obj = NULL) {
  ns <- getNamespace("dirichletprocess")
  
  if (is.null(dp_obj)) {
    # If no dp_obj provided, just check if C++ is available
    return(exists("_dirichletprocess_run_mcmc_cpp", where = ns))
  }

  if (!exists("_dirichletprocess_run_mcmc_cpp", where = ns)) {
    return(FALSE)
  }

  # Special case for mvnormal - needs specific functions
  if (inherits(dp_obj$mixingDistribution, "mvnormal")) {
    return(exists("conjugate_mvnormal_cluster_component_update_cpp", where = ns) &&
             exists("conjugate_mvnormal_cluster_parameter_update_cpp", where = ns))
  }

  # Supported types for unified MCMCRunner
  supported_types <- c("normal_inverse_gamma", "normal", "normalFixedVariance", "beta", "beta2",
                       "weibull", "exponential", "mvnormal", "mvnormal2")
  inherits(dp_obj$mixingDistribution, supported_types)
}

#' Run MCMC using C++ implementation
#' @param data Data matrix
#' @param mixing_dist_params Mixing distribution parameters
#' @param mcmc_params MCMC parameters
#' @return List with MCMC results
#' @keywords internal
run_mcmc_cpp <- function(data, mixing_dist_params, mcmc_params) {
  # Ensure required parameters
  if (!"m_auxiliary" %in% names(mcmc_params)) {
    mcmc_params$m_auxiliary <- 3  # Default for Algorithm 8
  }

  if (!"alpha" %in% names(mcmc_params)) {
    mcmc_params$alpha <- 1.0
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
#' @param dp_obj Dirichlet process object
#' @return List of parameters formatted for C++
#' @keywords internal
prepare_mixing_dist_params <- function(dp_obj) {
  md <- dp_obj$mixingDistribution

  if (inherits(md, "exponential")) {
    list(
      type = "exponential",
      alpha0 = md$priorParameters[1],
      beta0 = md$priorParameters[2]
    )
  } else if (inherits(md, "weibull")) {
    # Extract Weibull parameters
    list(
      type = "weibull",
      phi = md$priorParameters[1],      # Upper bound for alpha
      alpha0 = md$priorParameters[2],   # Shape for Gamma prior on 1/lambda
      beta0 = md$priorParameters[3],    # Rate for Gamma prior on 1/lambda
      hyper_a1 = ifelse(length(md$hyperPriorParameters) >= 1,
                        md$hyperPriorParameters[1], 6.0),
      hyper_a2 = ifelse(length(md$hyperPriorParameters) >= 2,
                        md$hyperPriorParameters[2], 2.0),
      hyper_b1 = ifelse(length(md$hyperPriorParameters) >= 3,
                        md$hyperPriorParameters[3], 1.0),
      hyper_b2 = ifelse(length(md$hyperPriorParameters) >= 4,
                        md$hyperPriorParameters[4], 0.5),
      mh_step_alpha = ifelse(!is.null(md$mhStepSize), md$mhStepSize[1], 0.1),
      mh_draws = ifelse(!is.null(dp_obj$mhDraws), dp_obj$mhDraws, 100)
    )
  } else if (inherits(md, "mvnormal")) {
    # Extract MVNormal parameters
    if (!is.null(md$priorParameters)) {
      pp <- md$priorParameters
      params <- list(
        type = "mvnormal",
        mu0 = as.numeric(pp$mu0),
        kappa0 = as.numeric(pp$kappa0),
        Lambda = as.matrix(pp$Lambda),
        nu = as.numeric(pp$nu)
      )
      
      # Add covariance model if specified
      if (!is.null(pp$covModel)) {
        params$covModel = as.character(pp$covModel)
      }
      
      return(params)
    } else {
      stop("MVNormal mixing distribution missing prior parameters")
    }
  } else if (inherits(md, "mvnormal2")) {
    # Extract MVNormal2 parameters
    if (!is.null(md$priorParameters)) {
      pp <- md$priorParameters
      list(
        type = "mvnormal2",
        mu0 = as.matrix(pp$mu0),
        sigma0 = as.matrix(pp$sigma0),
        phi0 = as.matrix(pp$phi0),
        nu0 = as.numeric(pp$nu0)
      )
    } else {
      stop("MVNormal2 mixing distribution missing prior parameters")
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
  } else if (inherits(md, "beta2")) {
    return(list(
      type = "beta2",
      gamma_prior = md$priorParameters[1],
      maxT = md$maxT,
      mh_step_size = md$mhStepSize,
      mh_draws = if (!is.null(md$mhDraws)) md$mhDraws else 250
    ))
  } else if (inherits(md, "normal_inverse_gamma") || inherits(md, "normal")) {
    # Gaussian parameters
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
  } else if (inherits(md, "normalFixedVariance")) {
    return(list(
      type = "normalFixedVariance",
      mu0 = md$priorParameters[1],
      sigma0 = md$priorParameters[2],
      sigma = md$sigma
    ))
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
#' @param dp_obj Dirichlet process object
#' @param its Number of iterations
#' @param updatePrior Whether to update prior parameters
#' @param n_burn Burn-in iterations
#' @param thin Thinning interval
#' @return List of MCMC parameters
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

#' Enable C++ implementations for specific samplers
#' @param enable Logical indicating whether to enable C++ samplers
#' @export
enable_cpp_samplers <- function(enable = TRUE) {
  if (missing(enable)) {
    # If no argument provided, return status (backward compatibility)
    return(invisible(exists("_dirichletprocess_run_mcmc_cpp", mode = "function")))
  }
  
  # Set option to force C++ usage
  options(dirichletprocess.force_cpp_samplers = enable)
  invisible(enable)
}

#' Enable C++ implementations for hierarchical models
#' @param enable Logical indicating whether to enable hierarchical C++ samplers
#' @export
enable_cpp_hierarchical_samplers <- function(enable = TRUE) {
  if (missing(enable)) {
    # If no argument provided, return status (backward compatibility)
    return(invisible(exists("_dirichletprocess_hierarchical_beta_fit_cpp", mode = "function")))
  }
  
  # Set option to force hierarchical C++ usage
  options(dirichletprocess.force_cpp_hierarchical = enable)
  invisible(enable)
}

#' Check if using C++ samplers
#' @export
using_cpp_samplers <- function() {
  ns <- getNamespace("dirichletprocess")
  
  # Check if forced via options
  force_cpp <- getOption("dirichletprocess.force_cpp_samplers", FALSE)
  if (force_cpp) {
    return(using_cpp() && exists("_dirichletprocess_run_mcmc_cpp", where = ns))
  }
  
  # Default behavior - check for existence of the compiled C++ function
  using_cpp() && exists("_dirichletprocess_run_mcmc_cpp", where = ns)
}

#' Check if using hierarchical C++ samplers
#' @export
using_cpp_hierarchical_samplers <- function() {
  ns <- getNamespace("dirichletprocess")
  
  # Check if forced via options
  force_hierarchical <- getOption("dirichletprocess.force_cpp_hierarchical", NULL)
  if (!is.null(force_hierarchical)) {
    if (force_hierarchical) {
      return(using_cpp() && exists("_dirichletprocess_hierarchical_beta_fit_cpp", where = ns))
    } else {
      return(FALSE)  # Force R implementation
    }
  }
  
  # Default behavior - C++ functions from Rcpp are stored as "list" mode, not "function"
  using_cpp() && exists("_dirichletprocess_hierarchical_beta_fit_cpp", where = ns)
}
