#' Update the component of the Dirichlet process
#'
#' Update the cluster assignment for each data point.
#'
#' @param dpObj Dirichlet Process object
#' @return Dirichlet process object with update components.
#'
#' @examples
#' dp <- DirichletProcessGaussian(rnorm(10))
#' dp <- ClusterComponentUpdate(dp)
#'
#' @export
ClusterComponentUpdate <- function(dpObj){
  UseMethod("ClusterComponentUpdate", dpObj)
}

#' @export
#' @rdname ClusterComponentUpdate
ClusterComponentUpdate.conjugate <- function(dpObj) {
  y <- dpObj$data
  n <- dpObj$n
  alpha <- dpObj$alpha

  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  mdObj <- dpObj$mixingDistribution

  pointsPerCluster <- dpObj$pointsPerCluster
  predictiveArray <- dpObj$predictiveArray

  for (i in seq_len(n)) {
    currentLabel <- clusterLabels[i]
    pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] - 1

    cluster_probs <- numeric(numLabels)

    for (j in 1:numLabels) {
      if (pointsPerCluster[j] > 0) {
        # Extract the parameters for cluster j, preserving dimensions
        single_cluster_params <- list()
        for (k in seq_along(clusterParams)) {
          param_dims <- dim(clusterParams[[k]])
          if (length(param_dims) == 3) {
            # For 3D arrays, extract the slice for cluster j
            single_cluster_params[[k]] <- array(
              clusterParams[[k]][, , j],
              dim = c(param_dims[1], param_dims[2], 1)
            )
          } else {
            # Fallback for other structures
            single_cluster_params[[k]] <- clusterParams[[k]][j]
          }
        }

        likelihood_val <- Likelihood(mdObj, y[i, , drop = FALSE], single_cluster_params)
        cluster_probs[j] <- pointsPerCluster[j] * as.numeric(likelihood_val[1])
      } else {
        cluster_probs[j] <- 0
      }
    }

    new_cluster_prob <- alpha * predictiveArray[i]
    probs <- c(cluster_probs, new_cluster_prob)

    probs[is.na(probs) | is.infinite(probs)] <- 0
    if (all(probs == 0)) {
      probs <- rep_len(1, length(probs))
    }

    newLabel <- sample.int(numLabels + 1, 1, prob = probs)

    dpObj$pointsPerCluster <- pointsPerCluster
    dpObj <- ClusterLabelChange(dpObj, i, newLabel, currentLabel)

    pointsPerCluster <- dpObj$pointsPerCluster
    clusterLabels <- dpObj$clusterLabels
    clusterParams <- dpObj$clusterParameters
    numLabels <- dpObj$numberClusters
  }

  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels
  return(dpObj)
}

#'@export
ClusterComponentUpdate.nonconjugate <- function(dpObj) {

  if (inherits(dpObj, "beta") && using_cpp_samplers()) {
    # Call the C++ implementation
    cpp_result <- nonconjugate_beta_cluster_component_update_cpp(dpObj)

    if (!is.null(cpp_result) && !isTRUE(cpp_result$stub_result)) {
      # If C++ implementation is complete and returns the updated dpObj structure
      dpObj$clusterLabels <- cpp_result$clusterLabels
      dpObj$pointsPerCluster <- cpp_result$pointsPerCluster
      dpObj$numberClusters <- cpp_result$numberClusters
      dpObj$clusterParameters <- cpp_result$clusterParameters
      return(dpObj)
    }
    # If cpp_result is NULL, it means C++ stub signaled to use R fallback (or an error occurred)
    # The warning from C++ will have already printed.
    # Proceed to R fallback logic below
  }

  # R fallback implementation (will run if C++ is disabled or C++ stub returns NULL)
  y <- dpObj$data
  n <- dpObj$n
  alpha <- dpObj$alpha

  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters

  mdObj <- dpObj$mixingDistribution
  m <- dpObj$m

  pointsPerCluster <- dpObj$pointsPerCluster

  aux <- vector("list", length(clusterParams))

  for (i in seq_len(n)) {
    currentLabel <- clusterLabels[i]

    pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] - 1

    # Determine the correct parameters for the cluster being emptied (or use prior if it was a singleton)
    current_params_for_empty_slot <- if(pointsPerCluster[currentLabel] == 0 && currentLabel <= dim(clusterParams[[1]])[3]) {
      list(
        mu = array(clusterParams[[1]][, , currentLabel], dim = c(1, 1, 1)),
        nu = array(clusterParams[[2]][, , currentLabel], dim = c(1, 1, 1))
      )
    } else {
      NULL
    }

    if (!is.null(current_params_for_empty_slot) && pointsPerCluster[currentLabel] == 0) {
      priorDraws_aux <- PriorDraw(mdObj, m - 1)
      aux[[1]] <- array(c(current_params_for_empty_slot[[1]], priorDraws_aux[[1]]), dim = c(1, 1, m))
      aux[[2]] <- array(c(current_params_for_empty_slot[[2]], priorDraws_aux[[2]]), dim = c(1, 1, m))
    } else {
      aux <- PriorDraw(mdObj, m)
    }

    cluster_probs <- numeric(numLabels)
    if (numLabels > 0) {
      for(k_idx in 1:numLabels) {
        if(pointsPerCluster[k_idx] > 0 && k_idx <= dim(clusterParams[[1]])[3]) {
          theta_k <- list()
          for (j in seq_along(clusterParams)) {
            param_dims <- dim(clusterParams[[j]])
            if (length(param_dims) == 3) {
              theta_k[[j]] <- array(clusterParams[[j]][,,k_idx], dim=c(param_dims[1], param_dims[2], 1))
            } else {
              theta_k[[j]] <- clusterParams[[j]][k_idx]
            }
          }
          cluster_probs[k_idx] <- pointsPerCluster[k_idx] * Likelihood(mdObj, y[i, , drop = FALSE], theta_k)
        } else {
          cluster_probs[k_idx] <- 0
        }
      }
    }

    aux_probs <- numeric(m)
    for(k_idx in 1:m) {
      theta_aux_k <- list()
      for (j in seq_along(aux)) {
        param_dims <- dim(aux[[j]])
        if (length(param_dims) == 3) {
          theta_aux_k[[j]] <- array(aux[[j]][,,k_idx], dim=c(param_dims[1], param_dims[2], 1))
        } else {
          theta_aux_k[[j]] <- aux[[j]][k_idx]
        }
      }
      aux_probs[k_idx] <- (alpha/m) * Likelihood(mdObj, y[i, , drop = FALSE], theta_aux_k)
    }

    probs <- c(cluster_probs, aux_probs)
    probs[is.na(probs) | !is.finite(probs)] <- 0

    if (all(probs == 0)) {
      probs <- rep_len(1, length(probs))
    }
    newLabel <- sample.int(length(probs), 1, prob = probs)

    # CRITICAL FIX: Update dpObj$pointsPerCluster BEFORE calling ClusterLabelChange
    dpObj$pointsPerCluster <- pointsPerCluster
    dpObj <- ClusterLabelChange(dpObj, i, newLabel, currentLabel, aux)

    # After a point is reassigned, the state of the clusters (number, labels, parameters)
    # might have changed. You MUST refresh the local variables from the returned dpObj
    # to ensure the next iteration of the loop has the most up-to-date information.
    pointsPerCluster <- dpObj$pointsPerCluster
    clusterLabels    <- dpObj$clusterLabels
    clusterParams    <- dpObj$clusterParameters
    numLabels        <- dpObj$numberClusters
  }

  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels
  return(dpObj)
}

#' @export
#' @rdname ClusterComponentUpdate
ClusterComponentUpdate.hierarchical <- function(dpObj){

  for(i in seq_along(dpObj$indDP)){
    dpObj$indDP[[i]] <- ClusterComponentUpdate(dpObj$indDP[[i]])
    dpObj$indDP[[i]] <- DuplicateClusterRemove(dpObj$indDP[[i]])
  }
  return(dpObj)
}
