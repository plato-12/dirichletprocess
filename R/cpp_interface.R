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

  # Call C++ function
  likelihood_cpp(mdObj, x, theta)
}

#' Register C++ Implementations
#'
#' This function registers the C++ implementations as alternatives to the R implementations.
#' It should be called during package loading.
#'
#' @keywords internal
register_cpp_implementations <- function() {
  # Create namespace environment for storing C++ implementations
  if (!exists("cpp_implementations", envir = parent.env(environment()))) {
    assign("cpp_implementations", new.env(), envir = parent.env(environment()))
  }

  # Register likelihood implementation
  cpp_env <- get("cpp_implementations", envir = parent.env(environment()))
  cpp_env$likelihood <- likelihood_cpp_wrapper
}

# Register implementations when package is loaded
.onLoad <- function(libname, pkgname) {
  register_cpp_implementations()
}
