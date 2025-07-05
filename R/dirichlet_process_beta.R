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
DirichletProcessBeta <- function(y, alphaPriors = c(2, 0.5),
                                 mhStepSize = c(0.1, 0.1), verbose = TRUE) {
  mdObj <- BetaMixtureCreate(priorParameters = c(2, 8),
                             mhStepSize = mhStepSize,
                             maxT = 1)

  dpObj <- DirichletProcessCreate(y, mdObj, alphaPriors)
  dpObj <- Initialise(dpObj, verbose = verbose)

  # Ensure cluster accounting is correct
  dpObj$pointsPerCluster <- as.numeric(table(factor(dpObj$clusterLabels,
                                                    levels = 1:dpObj$numberClusters)))

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
