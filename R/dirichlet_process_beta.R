#' Create a Dirichlet Process object with Beta mixing distribution
#'
#' @param y Data for which to be fitted
#' @param alphaPriors Alpha prior parameters for the DP concentration parameter
#' @param mhStepSize Metropolis-Hastings step size for parameter updates
#' @param verbose Logical indicating whether to print messages
#' @param cpp Logical. Use C++ implementation if TRUE, R implementation if FALSE. Default is FALSE.
#' @return Dirichlet process object with Beta mixing distribution
#' @export
DirichletProcessBeta <- function(y, alphaPriors = c(2, 0.5),
                                 mhStepSize = c(0.1, 0.1), verbose = TRUE,
                                 cpp = FALSE) {
  # Handle case where alphaPriors is a single value
  if (length(alphaPriors) == 1) {
    alphaPriors <- c(alphaPriors, 0.5)
  }
  mdObj <- BetaMixtureCreate(priorParameters = c(2, 8),
                             mhStepSize = mhStepSize,
                             maxT = 1)

  dpObj <- DirichletProcessCreate(y, mdObj, alphaPriors)
  dpObj <- Initialise(dpObj, verbose = verbose)

  # Ensure cluster accounting is correct
  dpObj$pointsPerCluster <- as.numeric(table(factor(dpObj$clusterLabels,
                                                    levels = 1:dpObj$numberClusters)))

  # Set cpp preference for this object
  if (cpp) {
    options(dirichletprocess.use_cpp = TRUE)
  } else {
    options(dirichletprocess.use_cpp = FALSE)
  }

  return(dpObj)
}

#' @export
#' @rdname Initialise
Initialise.beta <- function(dpObj, posterior = TRUE, m = 3, verbose = TRUE, numInitialClusters = 1, ...) {

  dpObj$m <- m
  dpObj$numberClusters <- 1
  dpObj$clusterLabels <- rep(1, dpObj$n)
  dpObj$pointsPerCluster <- c(dpObj$n)

  # Ensure parameters are properly structured as 3D arrays
  priorDraws <- PriorDraw(dpObj$mixingDistribution, 1)
  dpObj$clusterParameters <- list(
    mu = array(priorDraws$mu, dim = c(1, 1, 1)),
    nu = array(priorDraws$nu, dim = c(1, 1, 1))
  )

  dpObj$alpha <- dpObj$alphaPriorParameters[1] / dpObj$alphaPriorParameters[2]

  # Generate auxiliary parameters with proper structure
  dpObj$aux <- vector("list", m)
  for(i in seq_len(m)) {
    aux_draw <- PriorDraw(dpObj$mixingDistribution, 1)
    dpObj$aux[[i]] <- list(
      mu = array(aux_draw$mu, dim = c(1, 1, 1)),
      nu = array(aux_draw$nu, dim = c(1, 1, 1))
    )
  }

  if (verbose) {
    cat("Dirichlet process initialised.\n")
  }

  return(dpObj)
}
