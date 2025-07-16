#' The likelihood of the Dirichlet process object
#'
#' Calculate the likelihood of each data point with its parameter.
#'
#' @param dpobj The dirichletprocess object on which to calculate the likelihood.
#'
#' @export
LikelihoodDP <- function(dpobj){

  clusters_parameters <- dpobj$clusterParameters

  # For multivariate normal with pre-allocated arrays, we need to extract only active clusters
  if (inherits(dpobj, "mvnormal") && is.list(clusters_parameters)) {
    # Create a subset of parameters for only active clusters
    active_params <- list()
    for (i in seq_along(clusters_parameters)) {
      param_dims <- dim(clusters_parameters[[i]])
      if (length(param_dims) == 3 && param_dims[3] > dpobj$numberClusters) {
        # Extract only the active clusters
        if (dpobj$numberClusters == 1) {
          active_params[[i]] <- array(clusters_parameters[[i]][, , 1:dpobj$numberClusters, drop = FALSE],
                                      dim = c(param_dims[1], param_dims[2], dpobj$numberClusters))
        } else {
          active_params[[i]] <- clusters_parameters[[i]][, , 1:dpobj$numberClusters, drop = FALSE]
        }
      } else {
        active_params[[i]] <- clusters_parameters[[i]]
      }
    }
    clusters_parameters <- active_params
    names(clusters_parameters) <- names(dpobj$clusterParameters)
  }

  likelihoodValues <- vapply(seq_len(nrow(dpobj$data)),
                             function(i) {
                               lik <- Likelihood(dpobj$mixingDistribution,
                                                 dpobj$data[i, , drop=FALSE],
                                                 clusters_parameters)
                               # Ensure we only return values for active clusters
                               if (length(lik) > dpobj$numberClusters) {
                                 lik[1:dpobj$numberClusters]
                               } else {
                                 lik
                               }
                             },
                             numeric(dpobj$numberClusters))

  if (dpobj$numberClusters == 1) {
    likelihoodValues <- matrix(likelihoodValues, ncol = 1)
  } else {
    likelihoodValues <- t(likelihoodValues)
  }

  weight <- dpobj$pointsPerCluster / dpobj$n

  likelihoodValues <- likelihoodValues %*% weight

  return(likelihoodValues)
}
