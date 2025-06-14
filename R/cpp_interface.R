#' Interface to C++ Likelihood Implementation
#'
#' This function provides a bridge between the R S3 methods and the C++ implementation
#' of the likelihood functions.
#'
#' @param mdObj Mixing distribution object
#' @param x Data vector/matrix
#' @param theta List of parameters
#' @return Vector of likelihood values
#' @keywords internal
likelihood_cpp_wrapper <- function(mdObj, x, theta) {
  # Check if we have a distribution type
  if (!is.list(mdObj) || is.null(mdObj$distribution)) {
    return(UseMethod("Likelihood", mdObj))
  }

  dist_type <- mdObj$distribution

  # Handle exponential distribution directly
  if (dist_type == "exponential") {
    tryCatch({
      # Extract lambda parameter correctly
      if (is.list(theta) && length(theta) > 0) {
        lambda_array <- theta[[1]]
        if (is.array(lambda_array) || is.numeric(lambda_array)) {
          lambda <- as.numeric(lambda_array)[1]
          x_vec <- if (is.matrix(x)) as.vector(x) else as.numeric(x)
          return(exponential_likelihood_cpp(x_vec, lambda))
        }
      }
      # Fall back if parameter extraction fails
      return(UseMethod("Likelihood", mdObj))
    }, error = function(e) {
      return(UseMethod("Likelihood", mdObj))
    })
  }

  # Handle normal distribution
  if (dist_type == "normal") {
    tryCatch({
      if (is.list(theta) && length(theta) >= 2) {
        mu <- as.numeric(theta[[1]])[1]
        sigma <- as.numeric(theta[[2]])[1]
        x_vec <- if (is.matrix(x)) as.vector(x) else as.numeric(x)
        return(normal_likelihood_cpp(x_vec, mu, sigma))
      }
      return(UseMethod("Likelihood", mdObj))
    }, error = function(e) {
      return(UseMethod("Likelihood", mdObj))
    })
  }

  # Fall back to R implementation for other distributions
  return(UseMethod("Likelihood", mdObj))
}

#' Register C++ Implementations
#'
#' This function registers the C++ implementations as alternatives to the R implementations.
#' It should be called during package loading.
#'
#' @keywords internal
register_cpp_implementations <- function() {
  # Create namespace environment for storing C++ implementations
  pkg_env <- parent.env(environment())

  if (!exists("cpp_implementations", envir = pkg_env)) {
    assign("cpp_implementations", new.env(parent = pkg_env), envir = pkg_env)
  }

  # Register likelihood implementation when available
  cpp_env <- get("cpp_implementations", envir = pkg_env)
  cpp_env$likelihood <- likelihood_cpp_wrapper

  # Log that registration was attempted
  if (getOption("dirichletprocess.verbose", FALSE)) {
    message("C++ implementations registered")
  }
}

#' Get C++ Implementation Status
#'
#' Check which C++ implementations are available
#'
#' @return A named logical vector indicating which C++ implementations are available
#' @export
get_cpp_status <- function() {
  status <- list(
    likelihood = exists("likelihood_cpp"),
    normal_likelihood = exists("normal_likelihood_cpp"),
    exponential_likelihood = exists("exponential_likelihood_cpp"),  # Added
    mvnormal_likelihood = exists("mvnormal_likelihood_cpp"),
    benchmark_components = exists("benchmark_cpp_components"),
    memory_tracking = exists("get_memory_tracking")
  )
  return(status)
}

#' Check if C++ Implementation is Available
#'
#' @param func_name Name of the function to check
#' @return Logical indicating if the C++ implementation is available
#' @keywords internal
has_cpp_implementation <- function(func_name) {
  pkg_env <- parent.env(environment())

  if (!exists("cpp_implementations", envir = pkg_env)) {
    return(FALSE)
  }

  cpp_env <- get("cpp_implementations", envir = pkg_env)
  return(exists(func_name, envir = cpp_env))
}

# Register implementations when package is loaded
.onLoad <- function(libname, pkgname) {
  # Register C++ implementations
  register_cpp_implementations()

  # Set default options
  options(dirichletprocess.use_cpp = FALSE)
  options(dirichletprocess.verbose = FALSE)
}

.onAttach <- function(libname, pkgname) {
  # Check if C++ implementations are available
  cpp_status <- get_cpp_status()
  available_count <- sum(unlist(cpp_status))

  if (available_count > 0) {
    packageStartupMessage(
      sprintf("dirichletprocess: %d C++ implementations available", available_count)
    )

    if (getOption("dirichletprocess.verbose", FALSE)) {
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
  .Call("_dirichletprocess_run_mcmc_cpp",
        data = as.matrix(data),
        mixing_dist_params = mixing_dist_params,
        mcmc_params = mcmc_params,
        PACKAGE = "dirichletprocess")
}

#' Create mixing distribution parameters for C++
#' @keywords internal
prepare_mixing_dist_params <- function(dp_obj) {
  md <- dp_obj$mixing_distribution

  if (inherits(md, "normal_inverse_gamma")) {
    list(
      type = "gaussian",
      mu0 = md$priors$mu_0,
      kappa0 = md$priors$kappa_0,
      alpha0 = md$priors$alpha_0,
      beta0 = md$priors$beta_0
    )
  } else {
    stop("Mixing distribution not yet implemented in C++")
  }
}
