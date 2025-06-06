#' Initialise a Dirichlet process object
#'
#' Initialise a Dirichlet process object by assigning all the data points to a single cluster with a posterior or prior draw for parameters.
#'
#' @param dpObj A Dirichlet process object.
#' @param posterior TRUE/FALSE value for whether the cluster parameters should be from the posterior. If false then the values are from the prior.
#' @param m Number of auxiliary variables to use for a non-conjugate mixing distribution. Defaults to m=3. See \code{\link{ClusterComponentUpdate}} for more details on m.
#' @param verbose Logical flag indicating whether to output the acceptance ratio for non-conjugate mixtures.
#' @param numInitialClusters Number of clusters to initialise with.
#' @return A Dirichlet process object that has initial cluster allocations.
#' @export
Initialise <- function(dpObj, posterior = TRUE, m=3, verbose=TRUE, numInitialClusters = 1){
  UseMethod("Initialise", dpObj)
}

#' @export
Initialise.conjugate <- function(dpObj, posterior = TRUE, m=NULL, verbose=NULL, numInitialClusters = 1) {

  dpObj$clusterLabels <- rep_len(seq_len(numInitialClusters), length.out = dpObj$n)
  dpObj$numberClusters <- numInitialClusters
  dpObj$pointsPerCluster <- vapply(seq_len(numInitialClusters), function(x) sum(dpObj$clusterLabels == x), numeric(1))

  if (posterior && numInitialClusters == 1) {
    dpObj$clusterParameters <- PosteriorDraw(dpObj$mixingDistribution, dpObj$data, 1)
  } else {
    dpObj$clusterParameters <- PriorDraw(dpObj$mixingDistribution, numInitialClusters)
  }

  # For multivariate normal, ensure we have enough space for future clusters
  if (inherits(dpObj, "mvnormal")) {
    # Get current dimensions
    mu_dim <- dim(dpObj$clusterParameters$mu)
    sig_dim <- dim(dpObj$clusterParameters$sig)

    # Ensure we have at least enough slots for the data size or 50, whichever is larger
    min_slots <- max(50, dpObj$n, numInitialClusters * 10)

    if (mu_dim[3] < min_slots) {
      # Expand arrays
      d <- mu_dim[2]

      # Create new arrays with more space
      new_mu <- array(NA_real_, dim = c(1, d, min_slots))
      new_sig <- array(NA_real_, dim = c(d, d, min_slots))

      # Copy existing parameters
      new_mu[, , 1:mu_dim[3]] <- dpObj$clusterParameters$mu
      new_sig[, , 1:sig_dim[3]] <- dpObj$clusterParameters$sig

      # Fill remaining slots with prior draws
      if (mu_dim[3] < min_slots) {
        extra_params <- PriorDraw(dpObj$mixingDistribution, min_slots - mu_dim[3])
        new_mu[, , (mu_dim[3]+1):min_slots] <- extra_params$mu
        new_sig[, , (sig_dim[3]+1):min_slots] <- extra_params$sig
      }

      dpObj$clusterParameters$mu <- new_mu
      dpObj$clusterParameters$sig <- new_sig
    }
  }

  dpObj <- InitialisePredictive(dpObj)

  return(dpObj)
}

#'@export
Initialise.nonconjugate <- function(dpObj, posterior = TRUE, m = 3, verbose = TRUE, numInitialClusters=1) {

  # dpObj$clusterLabels <- 1:dpObj$n dpObj$numberClusters <- dpObj$n
  # dpObj$pointsPerCluster <- rep(1, dpObj$n) dpObj$clusterParameters <-
  # PosteriorDraw(dpObj$MixingDistribution, dpObj$data, dpObj$n)
  dpObj$clusterLabels <- rep(1, dpObj$n)
  dpObj$numberClusters <- 1
  dpObj$pointsPerCluster <- dpObj$n

  if (posterior) {
    post_draws <- PosteriorDraw(dpObj$mixingDistribution, dpObj$data, 1000)

    if (verbose)
      cat(paste("Accept Ratio: ",
                length(unique(c(post_draws[[1]])))/1000,
                "\n"))

    dpObj$clusterParameters <- lapply(post_draws, function(x) x[, , 1000, drop = FALSE])


    # dpObj$clusterParameters <- list(post_draws[[1]][, , 1000, drop = FALSE],
                                    # post_draws[[2]][, , 1000, drop = FALSE])
  } else {
    dpObj$clusterParameters <- PriorDraw(dpObj$mixingDistribution, 1)
  }

  dpObj$m <- m

  return(dpObj)
}


InitialisePredictive <- function(dpObj) UseMethod("InitialisePredictive", dpObj)

InitialisePredictive.conjugate <- function(dpObj) {

  dpObj$predictiveArray <- Predictive(dpObj$mixingDistribution, dpObj$data)

  return(dpObj)
}

InitialisePredictive.nonconjugate <- function(dpObj) {
  return(dpObj)
}



