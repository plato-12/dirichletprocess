#' Update the prior parameters of a mixing distribution
#'
#' @param mdObj Mixing Distribution Object
#' @param clusterParameters Current cluster parameters
#' @param n Number of samples
#' @return mdobj New Mixing Distribution object with updated cluster parameters
#' @export
PriorParametersUpdate <- function(mdObj, clusterParameters, n = 1){
  UseMethod("PriorParametersUpdate", mdObj)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.normal <- function(mdObj, clusterParameters, n = 1) {
  # For conjugate normal distributions, we typically don't update prior parameters
  # Return the original object unchanged
  if (getOption("dirichletprocess.verbose", FALSE)) {
    warning("Prior parameter update not implemented for conjugate normal distributions")
  }
  return(mdObj)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.conjugate <- function(mdObj, clusterParameters, n = 1) {
  # For conjugate distributions, we typically don't update prior parameters
  # Return the original object unchanged
  if (getOption("dirichletprocess.verbose", FALSE)) {
    warning("Prior parameter update not typically used for conjugate distributions")
  }
  return(mdObj)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.exponential <- function(mdObj, clusterParameters, n = 1) {
  # For exponential distributions, implement empirical Bayes update
  if (length(clusterParameters) == 0) {
    return(mdObj)
  }

  # Extract rate parameters from cluster parameters
  rates <- numeric(length(clusterParameters))
  for (i in seq_along(clusterParameters)) {
    if (is.list(clusterParameters[[i]]) && length(clusterParameters[[i]]) > 0) {
      rates[i] <- clusterParameters[[i]][[1]]
    } else if (is.numeric(clusterParameters[[i]])) {
      rates[i] <- clusterParameters[[i]][1]
    }
  }

  # Remove any invalid rates
  rates <- rates[rates > 0 & is.finite(rates)]

  if (length(rates) > 0) {
    # Update prior parameters based on observed rates
    mdObj$priorParameters[1] <- mean(rates)
    mdObj$priorParameters[2] <- var(rates) + mean(rates)^2
  }

  return(mdObj)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.weibull <- function(mdObj, clusterParameters, n = 1) {
  # For Weibull distributions, implement empirical Bayes update
  if (length(clusterParameters) == 0) {
    return(mdObj)
  }

  # Extract shape and scale parameters
  shapes <- numeric(length(clusterParameters))
  scales <- numeric(length(clusterParameters))

  for (i in seq_along(clusterParameters)) {
    if (is.list(clusterParameters[[i]]) && length(clusterParameters[[i]]) >= 2) {
      shapes[i] <- clusterParameters[[i]][[1]]
      scales[i] <- clusterParameters[[i]][[2]]
    }
  }

  # Remove invalid parameters
  valid_idx <- shapes > 0 & scales > 0 & is.finite(shapes) & is.finite(scales)
  shapes <- shapes[valid_idx]
  scales <- scales[valid_idx]

  if (length(shapes) > 0) {
    # Simple empirical Bayes update
    mdObj$priorParameters[1] <- mean(shapes)
    mdObj$priorParameters[2] <- mean(scales)
  }

  return(mdObj)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.beta <- function(mdObj, clusterParameters, n = 1) {
  # For beta distributions, implement empirical Bayes update
  if (length(clusterParameters) == 0) {
    return(mdObj)
  }

  # Extract alpha and beta parameters
  alphas <- numeric(length(clusterParameters))
  betas <- numeric(length(clusterParameters))

  for (i in seq_along(clusterParameters)) {
    if (is.list(clusterParameters[[i]]) && length(clusterParameters[[i]]) >= 2) {
      alphas[i] <- clusterParameters[[i]][[1]]
      betas[i] <- clusterParameters[[i]][[2]]
    }
  }

  # Remove invalid parameters
  valid_idx <- alphas > 0 & betas > 0 & is.finite(alphas) & is.finite(betas)
  alphas <- alphas[valid_idx]
  betas <- betas[valid_idx]

  if (length(alphas) > 0) {
    # Update prior parameters
    mdObj$priorParameters[1] <- mean(alphas)
    mdObj$priorParameters[2] <- mean(betas)
  }

  return(mdObj)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.default <- function(mdObj, clusterParameters, n = 1) {
  # Default implementation - return unchanged
  if (getOption("dirichletprocess.verbose", FALSE)) {
    warning("PriorParametersUpdate not implemented for this distribution type: ", class(mdObj))
  }
  return(mdObj)
}
