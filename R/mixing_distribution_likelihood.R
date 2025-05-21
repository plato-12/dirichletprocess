#' Mixing Distribution Likelihood
#'
#' Evaluate the Likelihood of some data \eqn{x} for some parameter \eqn{\theta}.
#'
#' @param mdObj Mixing Distribution
#' @param x Data
#' @param theta Parameters of distribution
#' @return Likelihood of the data
#' @export
Likelihood <- function(mdObj, x, theta) {
  if (using_cpp()) {
    tryCatch({
      return(likelihood_cpp_wrapper(mdObj, x, theta))
    }, error = function(e) {
      warning("C++ implementation failed, falling back to R: ", e$message)
    })
  }

  # Original implementation (falls back to this if C++ is not enabled or fails)
  UseMethod("Likelihood", mdObj)
}
