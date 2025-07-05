#' Draw from the posterior distribution
#'
#' @param mdObj Mixing Distribution
#' @param x Data
#' @param n Number of draws
#' @param ... For a non-conjugate distribution the starting parameters. Defaults to a draw from the prior distribution.
#' @return A sample from the posterior distribution
#' @export
PosteriorDraw <- function(mdObj, x, n = 1, ...){
  UseMethod("PosteriorDraw", mdObj)
}

#' @export
PosteriorDraw.nonconjugate <- function(mdObj, x, n = 1, ...) {

  if (missing(...) || is.null(list(...)$start_pos)) {
    # Try PenalisedLikelihood first, fall back to PriorDraw
    start_pos <- tryCatch({
      PenalisedLikelihood(mdObj, x)
    }, error = function(e) {
      PriorDraw(mdObj, 1)
    })
  } else {
    start_pos <- list(...)$start_pos
  }

  # Get MCMC samples
  mh_result <- MetropolisHastings(mdObj, x, start_pos, no_draws = n)

  # Extract parameter samples
  theta <- vector("list", length(mh_result$parameter_samples))

  for (i in seq_along(mh_result$parameter_samples)) {
    # Get the parameter array
    param_array <- mh_result$parameter_samples[[i]]

    # Ensure correct dimensions [1, 1, n]
    if (length(dim(param_array)) == 2) {
      # If 2D, reshape to 3D
      theta[[i]] <- array(param_array, dim = c(1, 1, n))
    } else if (length(dim(param_array)) == 3) {
      # Already 3D, just ensure it has the right number of samples
      theta[[i]] <- param_array[, , 1:n, drop = FALSE]
    } else {
      # Fallback: create array from values
      theta[[i]] <- array(as.numeric(param_array), dim = c(1, 1, n))
    }
  }

  return(theta)
}
