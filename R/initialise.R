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

    # Handle case where mu_dim might not have 3 dimensions (e.g., 1D data)
    if (is.null(mu_dim) || length(mu_dim) < 3) {
      # For 1D data, mu might be a vector or 2D array
      if (is.null(mu_dim)) {
        # It's a vector, convert to proper 3D array
        d <- 1
        n_clusters <- length(dpObj$clusterParameters$mu)
        dpObj$clusterParameters$mu <- array(dpObj$clusterParameters$mu, dim = c(1, d, n_clusters))
        dpObj$clusterParameters$sig <- array(dpObj$clusterParameters$sig, dim = c(d, d, n_clusters))
      } else if (length(mu_dim) == 2) {
        # It's a 2D array, add the third dimension
        d <- mu_dim[1]
        n_clusters <- mu_dim[2]
        dpObj$clusterParameters$mu <- array(dpObj$clusterParameters$mu, dim = c(1, d, n_clusters))
        # For constrained models, sig dimensions are different
        if (exists("priorParameters", dpObj$mixingDistribution) && 
            !is.null(dpObj$mixingDistribution$priorParameters$covModel) &&
            dpObj$mixingDistribution$priorParameters$covModel != "FULL") {
          # Keep sig as is for constrained models
        } else {
          dpObj$clusterParameters$sig <- array(dpObj$clusterParameters$sig, dim = c(d, d, n_clusters))
        }
      }
      # Update dimensions
      mu_dim <- dim(dpObj$clusterParameters$mu)
      sig_dim <- dim(dpObj$clusterParameters$sig)
    }

    # Ensure we have at least enough slots for the data size or 50, whichever is larger
    min_slots <- max(50, dpObj$n, numInitialClusters * 10)

    # Get current number of clusters from mu_dim, handling dimension issues
    current_clusters <- if (is.null(mu_dim) || length(mu_dim) < 3 || is.na(mu_dim[3])) {
      1  # Default to 1 cluster if dimension is problematic
    } else {
      mu_dim[3]
    }

    if (current_clusters < min_slots) {
      # Expand arrays
      d <- mu_dim[2]

      # Create new arrays with more space
      new_mu <- array(NA_real_, dim = c(1, d, min_slots))
      
      # For constrained models, sig dimensions are different
      if (exists("priorParameters", dpObj$mixingDistribution) && 
          !is.null(dpObj$mixingDistribution$priorParameters$covModel) &&
          dpObj$mixingDistribution$priorParameters$covModel != "FULL") {
        # Get number of parameters for this covariance model
        nParams <- dim(dpObj$clusterParameters$sig)[1]
        new_sig <- array(NA_real_, dim = c(nParams, min_slots))
        
        # Copy existing parameters
        new_mu[, , 1:current_clusters] <- dpObj$clusterParameters$mu
        new_sig[, 1:sig_dim[2]] <- dpObj$clusterParameters$sig
        
        # Fill remaining slots with prior draws
        if (current_clusters < min_slots) {
          extra_params <- PriorDraw(dpObj$mixingDistribution, min_slots - current_clusters)
          new_mu[, , (current_clusters+1):min_slots] <- extra_params$mu
          new_sig[, (sig_dim[2]+1):min_slots] <- extra_params$sig
        }
      } else {
        # Full covariance model
        new_sig <- array(NA_real_, dim = c(d, d, min_slots))
        
        # Copy existing parameters
        new_mu[, , 1:current_clusters] <- dpObj$clusterParameters$mu
        new_sig[, , 1:sig_dim[3]] <- dpObj$clusterParameters$sig
        
        # Fill remaining slots with prior draws
        if (current_clusters < min_slots) {
          extra_params <- PriorDraw(dpObj$mixingDistribution, min_slots - current_clusters)
          new_mu[, , (current_clusters+1):min_slots] <- extra_params$mu
          new_sig[, , (sig_dim[3]+1):min_slots] <- extra_params$sig
        }
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

# Covariance model-specific Initialise methods
#' @export
#' @rdname Initialise
Initialise.mvnormal.E <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  # Call base mvnormal initialise with covariance model handling
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.V <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.EII <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.VII <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.EEI <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.VEI <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.EVI <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.VVI <- function(dpObj, posterior = TRUE, m = NULL, verbose = NULL, numInitialClusters = 1) {
  return(Initialise.conjugate(dpObj, posterior, m, verbose, numInitialClusters))
}



