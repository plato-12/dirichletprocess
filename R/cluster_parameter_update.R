#'Update the cluster parameters of the Dirichlet process.
#'
#' Update the parameters of each individual cluster using all the data assigned to the particular cluster.
#' A sample is taken from the posterior distribution using a direct sample if the mixing distribution is conjugate or the Metropolis Hastings algorithm for non-conjugate mixtures.
#'
#'@param dpObj Dirichlet process object
#'@return Dirichlet process object with update cluster parameters
#'
#'@examples
#' dp <- DirichletProcessGaussian(rnorm(10))
#' dp <- ClusterParameterUpdate(dp)
#'
#'@export
ClusterParameterUpdate <- function(dpObj) UseMethod("ClusterParameterUpdate", dpObj)

#'@export
ClusterParameterUpdate.conjugate <- function(dpObj) {

  # Check for C++ implementation for MVNormal
  if (inherits(dpObj, "mvnormal") && using_cpp() &&
      exists("conjugate_mvnormal_cluster_parameter_update_cpp")) {
    return(ClusterParameterUpdate.mvnormal.cpp(dpObj))
  }

  y <- dpObj$data
  numLabels <- dpObj$numberClusters
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  mdobj <- dpObj$mixingDistribution

  # Check if numLabels is valid
  if (is.null(numLabels) || numLabels == 0) {
    return(dpObj)
  }

  for (i in 1:numLabels) {
    if (dpObj$pointsPerCluster[i] > 0) {
      pts <- y[which(clusterLabels == i), , drop = FALSE]
      post_draw <- PosteriorDraw(mdobj, pts)

      for (j in seq_along(clusterParams)) {
        param_dims <- dim(clusterParams[[j]])
        if (length(param_dims) == 3) {
          # FULL covariance model - 3D array
          clusterParams[[j]][, , i] <- post_draw[[j]]
        } else if (length(param_dims) == 2) {
          # Constrained covariance models - 2D array
          clusterParams[[j]][, i] <- post_draw[[j]]
        } else {
          # Single cluster case
          clusterParams[[j]][i] <- post_draw[[j]]
        }
      }
    }
  }

  dpObj$clusterParameters <- clusterParams
  return(dpObj)
}

#' @export
#' @rdname ClusterParameterUpdate
ClusterParameterUpdate.hierarchical <- function(dpObj) {
  # For hierarchical objects, update each individual DP
  for (i in seq_along(dpObj$indDP)) {
    dpObj$indDP[[i]] <- ClusterParameterUpdate(dpObj$indDP[[i]])
  }
  return(dpObj)
}

#' @export
#' @rdname ClusterParameterUpdate
ClusterParameterUpdate.nonconjugate <- function(dpObj) {

  if (inherits(dpObj, "beta") && using_cpp_samplers()) {
    cpp_result <- nonconjugate_beta_cluster_parameter_update_cpp(dpObj)
    if (!is.null(cpp_result)) {
      dpObj$clusterParameters <- cpp_result
      return(dpObj)
    }
  }

  for (i in seq_len(dpObj$numberClusters)) {
    cluster_data_indices <- dpObj$clusterLabels == i
    if (sum(cluster_data_indices) == 0) {
      next
    }
    cluster_data <- dpObj$data[cluster_data_indices, , drop = FALSE]

    # Prepare current parameters - handle different dimensions
    param1_dims <- dim(dpObj$clusterParameters[[1]])
    param2_dims <- dim(dpObj$clusterParameters[[2]])
    
    # Extract mu parameter
    if (length(param1_dims) == 3) {
      mu_param <- dpObj$clusterParameters[[1]][, , i]
    } else if (length(param1_dims) == 2) {
      mu_param <- dpObj$clusterParameters[[1]][, i]
    } else {
      mu_param <- dpObj$clusterParameters[[1]][i]
    }
    
    # Extract nu parameter  
    if (length(param2_dims) == 3) {
      nu_param <- dpObj$clusterParameters[[2]][, , i]
    } else if (length(param2_dims) == 2) {
      nu_param <- dpObj$clusterParameters[[2]][, i]
    } else {
      nu_param <- dpObj$clusterParameters[[2]][i]
    }
    
    current_params_list <- list(
      mu = array(mu_param, dim = c(1,1,1)),
      nu = array(nu_param, dim = c(1,1,1))
    )

    posterior_draw_samples <- PosteriorDraw(dpObj$mixingDistribution,
                                            cluster_data,
                                            n = dpObj$mhDraws,
                                            start_pos = current_params_list)

    # Handle different return formats from PosteriorDraw
    if (inherits(dpObj, "beta")) {
      # PosteriorDraw.beta returns list(mu=vector, nu=vector)
      # Extract the last sample from each
      mu_values <- posterior_draw_samples$mu
      nu_values <- posterior_draw_samples$nu

      # Take the last value from the MCMC chain - handle different dimensions
      if (length(param1_dims) == 3) {
        dpObj$clusterParameters[[1]][, , i] <- mu_values[length(mu_values)]
      } else if (length(param1_dims) == 2) {
        dpObj$clusterParameters[[1]][, i] <- mu_values[length(mu_values)]
      } else {
        dpObj$clusterParameters[[1]][i] <- mu_values[length(mu_values)]
      }
      
      if (length(param2_dims) == 3) {
        dpObj$clusterParameters[[2]][, , i] <- nu_values[length(nu_values)]
      } else if (length(param2_dims) == 2) {
        dpObj$clusterParameters[[2]][, i] <- nu_values[length(nu_values)]
      } else {
        dpObj$clusterParameters[[2]][i] <- nu_values[length(nu_values)]
      }
    } else {
      # Original logic for other distributions - handle different dimensions
      if (length(param1_dims) == 3) {
        dpObj$clusterParameters[[1]][, , i] <- posterior_draw_samples[[1]][,,dpObj$mhDraws, drop=FALSE]
      } else if (length(param1_dims) == 2) {
        dpObj$clusterParameters[[1]][, i] <- posterior_draw_samples[[1]][,dpObj$mhDraws, drop=FALSE]
      } else {
        dpObj$clusterParameters[[1]][i] <- posterior_draw_samples[[1]][dpObj$mhDraws]
      }
      
      if (length(param2_dims) == 3) {
        dpObj$clusterParameters[[2]][, , i] <- posterior_draw_samples[[2]][,,dpObj$mhDraws, drop=FALSE]
      } else if (length(param2_dims) == 2) {
        dpObj$clusterParameters[[2]][, i] <- posterior_draw_samples[[2]][,dpObj$mhDraws, drop=FALSE]
      } else {
        dpObj$clusterParameters[[2]][i] <- posterior_draw_samples[[2]][dpObj$mhDraws]
      }
    }
  }
  return(dpObj)
}
