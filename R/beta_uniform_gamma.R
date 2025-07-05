#' Create a Beta mixing distribution.
#'
#' See \code{\link{DirichletProcessBeta}} for the default prior and hyper prior distributions.
#'
#' @param priorParameters The prior parameters for the base measure.
#' @param mhStepSize The Metropolis Hastings step size. A numeric vector of length 2.
#' @param maxT The upper bound of the Beta distribution. Defaults to 1 for the standard Beta distribution.
#' @param hyperPriorParameters The parameters for the hyper prior.
#' @return A mixing distribution object.
#' @export
BetaMixtureCreate <- function(priorParameters = c(2, 8), mhStepSize = c(1, 1), maxT = 1,
                              hyperPriorParameters = c(1, 0.125)) {

  mdObj <- MixingDistribution("beta",
                              priorParameters, "nonconjugate",
                              mhStepSize, hyperPriorParameters)
  mdObj$maxT <- maxT
  return(mdObj)
}

#' @export
#' @rdname Likelihood
Likelihood.beta <- function(mdObj, x, theta) {
  maxT <- mdObj$maxT
  x <- as.vector(x, "numeric")
  mu <- as.numeric(theta[[1]][, , , drop = TRUE])
  tau <- as.numeric(theta[[2]][, , , drop = TRUE])

  # Ensure we have values
  if (length(mu) == 0 || length(tau) == 0) {
    return(numeric(length(x)))
  }

  # Recycle parameters if needed
  mu <- rep_len(mu, length(x))
  tau <- rep_len(tau, length(x))

  # Calculate likelihood
  y <- numeric(length(x))
  for (i in seq_along(x)) {
    # Validate parameters
    if (is.na(mu[i]) || is.na(tau[i]) || mu[i] <= 0 || mu[i] >= maxT || tau[i] <= 0) {
      y[i] <- 1e-300
      next
    }

    a <- (mu[i] * tau[i]) / maxT
    b <- (1 - mu[i]/maxT) * tau[i]

    # Ensure valid beta parameters
    if (a <= 0 || b <= 0 || !is.finite(a) || !is.finite(b)) {
      y[i] <- 1e-300
      next
    }

    # Calculate likelihood
    if (x[i] >= 0 && x[i] <= maxT) {
      y[i] <- (1/maxT) * dbeta(x[i]/maxT, a, b)
    } else {
      y[i] <- 1e-300
    }
  }

  return(as.numeric(y))
}

#' @export
#' @rdname PriorDraw
PriorDraw.beta <- function(mdObj, n = 1) {

  priorParameters <- mdObj$priorParameters
  mu <- runif(n, 0, mdObj$maxT)
  nu <- 1/rgamma(n, shape = priorParameters[1], rate = priorParameters[2])

  theta <- list(mu = array(mu, c(1, 1, n)), nu = array(nu, c(1, 1, n)))
  return(theta)
}

#' @export
#' @rdname PriorDensity
PriorDensity.beta <- function(mdObj, theta) {

  priorParameters <- mdObj$priorParameters
  mu <- theta[[1]]
  nu <- theta[[2]]

  muDensity <- dunif(mu, 0, mdObj$maxT)

  nuDensity <- dgamma(1/nu, priorParameters[1], priorParameters[2]) * (1/nu^2)

  if(is.infinite(nuDensity) | is.na(nuDensity)){
    nuDensity <- 1e-10 # Return a very small number instead of Inf or NA
  }

  thetaDensity <- muDensity * nuDensity
  return(as.numeric(thetaDensity))
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.beta <- function(mdObj, clusterParameters, n = 1) {

  hyperPriorParameters <- mdObj$hyperPriorParameters
  priorParameters <- mdObj$priorParameters

  numClusters <- dim(clusterParameters[[1]])[3]

  posteriorShape <- hyperPriorParameters[1] + priorParameters[1] * numClusters
  posteriorRate <- hyperPriorParameters[2] + sum(1/clusterParameters[[2]])

  newGamma <- rgamma(n, posteriorShape, posteriorRate)

  newPriorParameters <- matrix(c(priorParameters[1], newGamma), ncol = 2)
  mdObj$priorParameters <- newPriorParameters

  return(mdObj)
}

#' @export
#' @rdname MhParameterProposal
MhParameterProposal.beta <- function(mdObj, old_params) {

  mhStepSize <- mdObj$mhStepSize

  new_params <- old_params

  # Extract current values
  old_mu <- as.numeric(old_params[[1]])
  old_nu <- as.numeric(old_params[[2]])

  # Propose new mu
  new_mu <- old_mu + mhStepSize[1] * rnorm(1, 0, 2.4)
  if (new_mu > mdObj$maxT || new_mu < 0) {
    new_mu <- old_mu
  }

  # Propose new nu (ensure positive)
  new_nu <- abs(old_nu + mhStepSize[2] * rnorm(1, 0, 2.4))

  # Return in proper format
  new_params[[1]] <- array(new_mu, dim = c(1, 1, 1))
  new_params[[2]] <- array(new_nu, dim = c(1, 1, 1))

  return(new_params)
}

#' @export
#' @rdname PenalisedLikelihood
PenalisedLikelihood.beta <- function(mdObj, x) {
  if (length(x) == 0) {
    return(PriorDraw(mdObj, 1))
  }

  x <- as.numeric(x)
  x <- x[x > 0 & x < mdObj$maxT]  # Remove boundary values

  if (length(x) == 0) {
    return(PriorDraw(mdObj, 1))
  }

  # Method of moments estimation
  x_norm <- x / mdObj$maxT
  x_mean <- mean(x_norm)
  x_var <- var(x_norm)

  # Handle edge cases
  if (is.na(x_var) || x_var < 1e-10) {
    x_var <- 0.01
  }

  if (x_mean <= 0.01) x_mean <- 0.01
  if (x_mean >= 0.99) x_mean <- 0.99

  # Calculate parameters
  common <- x_mean * (1 - x_mean) / x_var - 1
  if (common <= 0) {
    # Fallback to prior
    return(PriorDraw(mdObj, 1))
  }

  mu_est <- x_mean * mdObj$maxT
  tau_est <- common

  # Return in the expected format
  return(list(
    mu = array(mu_est, dim = c(1, 1, 1)),
    nu = array(tau_est, dim = c(1, 1, 1))
  ))
}
