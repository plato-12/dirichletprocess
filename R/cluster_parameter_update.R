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

  y <- dpObj$data
  numLabels <- dpObj$numberClusters
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  mdobj <- dpObj$mixingDistribution

  for (i in 1:numLabels) {
    # Safeguard: only update clusters with points
    if (dpObj$pointsPerCluster[i] > 0) {
      pts <- y[which(clusterLabels == i), , drop = FALSE]
      post_draw <- PosteriorDraw(mdobj, pts)

      for (j in seq_along(clusterParams)) {
        clusterParams[[j]][, , i] <- post_draw[[j]]
      }
    }
  }

  dpObj$clusterParameters <- clusterParams
  return(dpObj)
}

#' @export
#' @rdname ClusterParameterUpdate
ClusterParameterUpdate.nonconjugate <- function(dpObj) {

  if (inherits(dpObj, "beta") && using_cpp_samplers()) {
    # Call the C++ implementation STUB
    cpp_result <- nonconjugate_beta_cluster_parameter_update_cpp(dpObj)

    if (!is.null(cpp_result)) {
      # If C++ were fully implemented, it would return the updated dpObj
      # For a stub returning R_NilValue (which becomes NULL in R), this block is skipped
      # dpObj <- cpp_result # Or assign components like in component_update
      # return(dpObj)

      # If the stub just returns dp_list, this line tries to re-assign
      # This was likely the source of "length 11" if cpp_result was not dpObj
      # but some other list created by Rcpp by mistake.
      # Since the stub now returns R_NilValue, this will be skipped.
      dpObj$clusterParameters <- cpp_result$clusterParameters
      return(dpObj)
    }
    # If cpp_result is NULL, proceed to R fallback
  }

  # R fallback implementation
  for (i in seq_len(dpObj$numberClusters)) {
    cluster_data_indices <- dpObj$clusterLabels == i
    if (sum(cluster_data_indices) == 0) { # No data points in this cluster
      next
    }
    cluster_data <- dpObj$data[cluster_data_indices, , drop = FALSE]

    # Current parameters for cluster i (to be used as start for MH)
    current_params_list <- list(
      mu = array(dpObj$clusterParameters[[1]][, , i], dim = c(1,1,1)),
      nu = array(dpObj$clusterParameters[[2]][, , i], dim = c(1,1,1))
    )

    # PosteriorDraw for non-conjugate returns a list of samples (mu, nu)
    # We need to take the last sample as the updated parameter
    posterior_draw_samples <- PosteriorDraw(dpObj$mixingDistribution,
                                            cluster_data,
                                            current_params_list, # Start for MH
                                            num_draws = dpObj$mhDraws) # Assuming PosteriorDraw for non-conj takes num_draws

    # Update the parameters for cluster i with the last sample from MH
    dpObj$clusterParameters[[1]][, , i] <- posterior_draw_samples$mu[,,dpObj$mhDraws, drop=FALSE]
    dpObj$clusterParameters[[2]][, , i] <- posterior_draw_samples$nu[,,dpObj$mhDraws, drop=FALSE]
  }
  return(dpObj)
}

cluster_parameter_update <- function(mdobj, data, clusters, params){

  uniqueClusters <- unique(clusters)

  newParams <- lapply(uniqueClusters, function(i){
    updateData <- data[clusters==i, ,drop=F]
    newParam <- PosteriorDraw(mdobj, updateData)
    return(newParam)

  } )

  return(newParams)
}
