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
  # Convert data to appropriate format
  if (is.matrix(x)) {
    x <- as.vector(x)
  }

  # For now, we'll fall back to the R implementation
  # This will be replaced with actual C++ calls once they're implemented
  tryCatch({
    # Try to call C++ implementation if available
    if (exists("likelihood_cpp")) {
      return(likelihood_cpp(mdObj, x, theta))
    } else {
      # Fall back to R implementation
      return(UseMethod("Likelihood", mdObj))
    }
  }, error = function(e) {
    # If C++ implementation fails, fall back to R
    return(UseMethod("Likelihood", mdObj))
  })
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
