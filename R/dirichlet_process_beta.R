#' Create a Dirichlet Process object with Beta mixing distribution
#'
#' @param y Data for which to be fitted
#' @param maxY Upper bound for the Beta distribution
#' @param g0Priors Prior parameters for the base measure (alpha, beta for the prior)
#' @param alphaPriors Alpha prior parameters for the DP concentration parameter
#' @param mhStepSize Metropolis-Hastings step size for parameter updates
#' @param hyperPriorParameters Hyper prior parameters
#' @param verbose Logical indicating whether to print messages
#' @return Dirichlet process object with Beta mixing distribution
#' @export
DirichletProcessBeta <- function(y,
                                 maxY = NULL,
                                 g0Priors = c(2, 8),
                                 alphaPriors = c(2, 4),
                                 mhStepSize = c(1, 1),
                                 hyperPriorParameters = c(1, 0.125),
                                 verbose = TRUE) {

  if (is.null(maxY)) {
    maxY <- max(y) + 0.01  # Small buffer to ensure all data is included
  }

  mdObj <- BetaMixtureCreate(
    priorParameters = g0Priors,
    mhStepSize = mhStepSize,
    maxT = maxY,
    hyperPriorParameters = hyperPriorParameters
  )

  dpObj <- DirichletProcessCreate(
    y,
    mdObj,
    alphaPriors = alphaPriors,
    verbose = verbose
  )

  dpObj <- Initialise(dpObj)

  return(dpObj)
}
