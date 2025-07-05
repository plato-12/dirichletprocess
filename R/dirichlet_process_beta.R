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
Initialise.beta <- function(dpObj, posterior = TRUE, verbose = TRUE, ...) {

  # Ensure all points start in cluster 1
  dpObj$clusterLabels <- rep(1, dpObj$n)
  dpObj$numberClusters <- 1
  dpObj$pointsPerCluster <- numeric(dpObj$n)
  dpObj$pointsPerCluster[1] <- dpObj$n

  # Initialize parameters
  if (posterior) {
    cluster_data <- matrix(dpObj$data, ncol = 1)
    post_draws <- PosteriorDraw(dpObj$mixingDistribution, cluster_data, n = 1)
    dpObj$clusterParameters <- list(
      mu = as.numeric(post_draws$mu),
      nu = as.numeric(post_draws$nu)
    )
  } else {
    prior_draws <- PriorDraw(dpObj$mixingDistribution, 1)
    dpObj$clusterParameters <- list(
      mu = as.numeric(prior_draws$mu),
      nu = as.numeric(prior_draws$nu)
    )
  }

  # Initialize auxiliary parameters for non-conjugate
  dpObj$m <- 3
  dpObj$aux <- vector("list", dpObj$m)
  for (j in seq_len(dpObj$m)) {
    aux_params <- PriorDraw(dpObj$mixingDistribution, 1)
    dpObj$aux[[j]] <- list(
      mu = as.numeric(aux_params$mu),
      nu = as.numeric(aux_params$nu)
    )
  }

  if (verbose) {
    cat("Initialised Dirichlet process with 1 cluster\n")
  }

  return(dpObj)
}
