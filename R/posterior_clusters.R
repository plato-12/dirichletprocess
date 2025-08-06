#' Generate the posterior clusters of a Dirichlet Process
#'
#' Using the stick breaking representation the user can draw the posterior clusters and weights for a fitted Dirichlet Process.
#' See also \code{\link{PosteriorFunction}}.
#'
#' @param dpobj Fitted Dirichlet process
#' @param ind Index for which the posterior will be drawn from. Defaults to the last iteration of the fit.
#' @return A list with the weights and cluster parameters that form the posterior of the Dirichlet process.
#'
#' @examples
#' y <- rnorm(10)
#' dp <- DirichletProcessGaussian(y)
#' dp <- Fit(dp, 5)
#' postClusters <- PosteriorClusters(dp)
#'
#' @export
PosteriorClusters <- function(dpobj, ind) UseMethod("PosteriorClusters", dpobj)

#' @export
PosteriorClusters.dirichletprocess <- function(dpobj, ind) {

  if (!missing(ind)) {
    pointsPerCluster <- dpobj$weightsChain[[ind]] * dpobj$n
    alpha <- dpobj$alphaChain[ind]
    clusterParams <- dpobj$clusterParametersChain[[ind]]
  } else {
    pointsPerCluster <- dpobj$pointsPerCluster
    alpha <- dpobj$alpha
    clusterParams <- dpobj$clusterParameters
  }

  numLabels <- length(pointsPerCluster)
  mdobj <- dpobj$mixingDistribution

  # Remove zero clusters to avoid invalid arguments in rdirichlet
  non_zero_clusters <- pointsPerCluster > 0
  active_pointsPerCluster <- pointsPerCluster[non_zero_clusters]
  
  # If no active clusters, create a minimum viable cluster
  if (length(active_pointsPerCluster) == 0) {
    active_pointsPerCluster <- c(1)
  }
  
  # Ensure alpha is numeric
  if (!is.numeric(alpha)) {
    alpha <- as.numeric(alpha)
  }
  
  dirichlet_draws <- gtools::rdirichlet(1, c(active_pointsPerCluster, alpha))
  numBreaks <- ceiling(alpha + numLabels) * 20 + 5

  sticks <- StickBreaking(alpha + numLabels, numBreaks)
  active_numLabels <- length(active_pointsPerCluster)
  sticks <- sticks * dirichlet_draws[active_numLabels + 1]

  # Build the full sticks vector including zeros for empty clusters
  full_sticks <- numeric(numLabels)
  full_sticks[non_zero_clusters] <- dirichlet_draws[1:active_numLabels]
  sticks <- c(full_sticks, sticks)
  # postParams <- rbind(clusterParams, PriorDraw(mdobj, numBreaks))

  #n_smps <- numBreaks + numLabels

  PriorDraws <- PriorDraw(mdobj, numBreaks)
  postParams <- list()

  # For normal distributions, clusterParams contains the actual parameter data
  # and we need to construct synthetic parameter arrays for stick-breaking
  if (inherits(mdobj, "normal") && inherits(mdobj, "conjugate")) {
    # For normal conjugate case, create synthetic parameter structure
    # that matches what the plotting functions expect
    
    # Get the parameter names from PriorDraws
    param_names <- names(PriorDraws)
    
    if (!is.null(param_names) && length(param_names) == 2) {
      # For normal distribution: mu and sigma parameters
      # Create arrays that combine existing cluster data with prior draws
      for (i in seq_along(param_names)) {
        param_name <- param_names[i]
        
        # PriorDraws has the right structure: [1, 1, numBreaks]
        # We need to create a compatible structure for existing clusters
        cluster_array <- array(0, dim = c(1, 1, numLabels))
        
        # Fill the cluster array with actual parameter values (if available)
        # For normal distribution, we use synthetic values based on data
        if (numLabels > 0) {
          if (param_name == "mu") {
            # Use cluster means as synthetic mu values
            for (k in seq_len(numLabels)) {
              cluster_array[1, 1, k] <- mean(dpobj$data[dpobj$clusterLabels == k])
            }
          } else { # sigma
            # Use cluster standard deviations as synthetic sigma values
            for (k in seq_len(numLabels)) {
              cluster_data <- dpobj$data[dpobj$clusterLabels == k]
              cluster_array[1, 1, k] <- if(length(cluster_data) > 1) sd(cluster_data) else 1.0
            }
          }
        }
        
        # Combine cluster parameters with prior draws
        postParams[[i]] <- array(c(cluster_array, PriorDraws[[param_name]]),
                                 dim = c(1, 1, numBreaks + numLabels))
      }
      names(postParams) <- param_names
    }
  } else {
    # For other distributions, use the original logic with better error handling
    param_names <- names(clusterParams)
    if (is.null(param_names)) {
      param_names <- names(PriorDraws)
    }
    
    # Ensure we have valid names and matching structure
    if (is.null(param_names) || length(param_names) != length(clusterParams)) {
      # Fallback: use numeric indices but check bounds
      for (i in seq_along(clusterParams)) {
        if (i <= length(PriorDraws)) {
          postParams[[i]] <- array(c(clusterParams[[i]], PriorDraws[[i]]),
                                   dim = c(dim(PriorDraws[[i]])[1:2],
                                           numBreaks + numLabels))
        } else {
          # If PriorDraws is shorter, just use cluster params
          postParams[[i]] <- clusterParams[[i]]
        }
      }
    } else {
      # Use names to match parameters correctly
      for (i in seq_along(param_names)) {
        param_name <- param_names[i]
        if (param_name %in% names(PriorDraws)) {
          postParams[[i]] <- array(c(clusterParams[[i]], PriorDraws[[param_name]]),
                                   dim = c(dim(PriorDraws[[param_name]])[1:2],
                                           numBreaks + numLabels))
        } else {
          # If parameter not found in PriorDraws, just use cluster params
          postParams[[i]] <- clusterParams[[i]]
        }
      }
      names(postParams) <- param_names
    }
  }


  # smps <- sample.int(n_smps, replace = T, prob = sticks)
  # smpTable <- data.frame((table(smps)/n_smps))
  #   retParams <- list()
  # for (i in seq_along(clusterParams)) {
  #   retParams[[i]] <- postParams[[i]][, , smpTable$smp, drop = FALSE]
  # }
  # returnList <- list(weights = smpTable$Freq, params = retParams)
  #

  returnList <- list(weights=sticks, params=postParams)

  return(returnList)
}



