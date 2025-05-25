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

  # Check if numLabels is valid
  if (is.null(numLabels) || numLabels == 0) {
    return(dpObj)
  }

  for (i in 1:numLabels) {
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

    current_params_list <- list(
      mu = array(dpObj$clusterParameters[[1]][, , i], dim = c(1,1,1)),
      nu = array(dpObj$clusterParameters[[2]][, , i], dim = c(1,1,1))
    )

    posterior_draw_samples <- PosteriorDraw(dpObj$mixingDistribution,
                                            cluster_data,
                                            n = dpObj$mhDraws,
                                            start_pos = current_params_list)

    dpObj$clusterParameters[[1]][, , i] <- posterior_draw_samples[[1]][,,dpObj$mhDraws, drop=FALSE]
    dpObj$clusterParameters[[2]][, , i] <- posterior_draw_samples[[2]][,,dpObj$mhDraws, drop=FALSE]
  }
  return(dpObj)
}
