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
  mu <- theta[[1]][, , , drop = TRUE]
  tau <- theta[[2]][, , , drop = TRUE]


  a <- (mu * tau)/maxT
  b <- (1 - mu/maxT) * tau

  y <- 1/maxT * dbeta(x/maxT, a, b)

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
PenalisedLikelihood.beta <- function(mdObj, x){

  optimStartParams <- c(mdObj$maxT/2, 2)

  optimParams <- tryCatch(optim(optimStartParams, function(params){

    # Ensure params are in a valid range for log calculation
    params_mu <- params[1]
    params_nu <- params[2]
    if(params_mu <= 0 || params_mu >= mdObj$maxT || params_nu <= 0) return(1e30)


    ll <- sum(log(Likelihood(mdObj, x, VectorToArray(params))))
    # Ensure PriorDensity does not return 0 or negative before taking log
    prior_dens <- PriorDensity(mdObj, VectorToArray(params))
    if(prior_dens <= 0) return(1e30)
    ll <- ll + log(prior_dens)


    if (is.infinite(ll) || is.na(ll)) ll <- -1e30 # Use a large negative finite number

    return(-ll)
  }, method="L-BFGS-B", lower=c(1e-6,1e-6), upper=c(mdObj$maxT - 1e-6, Inf)), error = function(e) list(par=optimStartParams))


  optimParamsRet <- VectorToArray(optimParams$par)

  return(optimParamsRet)
}
