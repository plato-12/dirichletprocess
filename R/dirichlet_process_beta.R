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

  # FIXED: Pass alphaPriorParameters instead of alphaPriors
  # FIXED: Remove verbose parameter as DirichletProcessCreate doesn't accept it
  dpObj <- DirichletProcessCreate(
    y,
    mdObj,
    alphaPriorParameters = alphaPriors  # Changed from alphaPriors to alphaPriorParameters
  )

  # Store verbose for later use if needed
  dpObj$verbose <- verbose

  dpObj <- Initialise(dpObj, verbose = verbose)

  return(dpObj)
}

#' @export
Initialise.beta <- function(mdObj, dpObj, posterior = TRUE, verbose = TRUE, ...) {

  dpObj <- NextMethod()

  # Ensure proper initialization of cluster components
  if (dpObj$n > 0) {
    # Initialize all points to cluster 1
    dpObj$clusterLabels <- rep(1, dpObj$n)
    dpObj$pointsPerCluster <- dpObj$n
    dpObj$numberClusters <- 1

    # Initialize cluster parameters
    if (posterior) {
      dpObj$clusterParameters <- PosteriorDraw(dpObj$mixingDistribution,
                                               dpObj$data,
                                               dpObj$numberClusters)
    } else {
      dpObj$clusterParameters <- PriorDraw(dpObj$mixingDistribution,
                                           dpObj$numberClusters)
    }
  }

  if (verbose){
    cat("Initialised Dirichlet process with a", dpObj$mixingDistribution$distribution,
        "mixing distribution.\n")
  }
  return(dpObj)
}
