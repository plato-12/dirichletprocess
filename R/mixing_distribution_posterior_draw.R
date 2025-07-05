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

#' @export
PosteriorDraw.beta <- function(mdObj, x, n = 1, ...) {

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

  # Extract and return with proper names
  if (length(mh_result$parameter_samples) >= 2) {
    mu_samples <- mh_result$parameter_samples[[1]]
    nu_samples <- mh_result$parameter_samples[[2]]

    # Ensure correct dimensions
    if (length(dim(mu_samples)) < 3) {
      mu_samples <- array(mu_samples, dim = c(1, 1, n))
    }
    if (length(dim(nu_samples)) < 3) {
      nu_samples <- array(nu_samples, dim = c(1, 1, n))
    }

    return(list(mu = mu_samples, nu = nu_samples))
  } else {
    # Fallback to empty arrays
    return(list(
      mu = array(numeric(0), dim = c(1, 1, 0)),
      nu = array(numeric(0), dim = c(1, 1, 0))
    ))
  }
}
