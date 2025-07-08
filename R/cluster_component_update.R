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

  # Check for C++ implementation for MVNormal
  if (inherits(dpObj, "mvnormal") && using_cpp() &&
      exists("conjugate_mvnormal_cluster_component_update_cpp")) {
    return(ClusterComponentUpdate.mvnormal.cpp(dpObj))
  }

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
  }

  # R implementation fallback
  y <- dpObj$data
  n <- nrow(y)
  alpha <- dpObj$alpha
  m <- dpObj$m

  pointsPerCluster <- dpObj$pointsPerCluster
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  mdObj <- dpObj$mixingDistribution

  # Algorithm 8 from Neal (2000)
  for (i in seq_len(n)) {
    currentLabel <- clusterLabels[i]

    # Temporarily remove the point from its current cluster
    pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] - 1

    if (pointsPerCluster[currentLabel] == 0) {
      # If cluster is now empty, remove it
      numLabels <- numLabels - 1

      if (numLabels == 0) {
        # Edge case: all points were in one cluster, now empty
        # We'll reassign this point to a new cluster (from auxiliary)
        pointsPerCluster <- numeric(0)
        clusterParams <- lapply(clusterParams, function(x) array(numeric(0), dim = c(dim(x)[1], dim(x)[2], 0)))
      } else {
        # Remove the empty cluster
        pointsPerCluster <- pointsPerCluster[-currentLabel]

        # Remove parameters for empty cluster
        for (j in seq_along(clusterParams)) {
          clusterParams[[j]] <- clusterParams[[j]][, , -currentLabel, drop = FALSE]
        }

        # Adjust labels for all points
        for (k in seq_len(n)) {
          if (clusterLabels[k] > currentLabel) {
            clusterLabels[k] <- clusterLabels[k] - 1
          }
        }
      }
    }

    # Draw m auxiliary parameters from prior
    aux <- PriorDraw(mdObj, m)

    # Calculate probabilities for existing clusters
    if (numLabels > 0) {
      cluster_probs <- numeric(numLabels)
      for (j in seq_len(numLabels)) {
        cluster_count <- pointsPerCluster[j]
        cluster_data_idx <- which(clusterLabels == j)
        if (length(cluster_data_idx) > 0) {
          tempThetaJ <- list()
          for (k in seq_along(clusterParams)) {
            tempThetaJ[[k]] <- array(clusterParams[[k]][, , j],
                                     dim = c(dim(clusterParams[[k]])[1],
                                             dim(clusterParams[[k]])[2], 1))
          }
          cluster_probs[j] <- cluster_count * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
        } else {
          cluster_probs[j] <- 0
        }
      }
    } else {
      cluster_probs <- numeric(0)
    }

    # Calculate probabilities for auxiliary clusters
    aux_probs <- numeric(m)
    for (j in seq_len(m)) {
      tempThetaJ <- list()
      for (k in seq_along(aux)) {
        tempThetaJ[[k]] <- array(aux[[k]][, , j],
                                 dim = c(dim(aux[[k]])[1],
                                         dim(aux[[k]])[2], 1))
      }
      aux_probs[j] <- (alpha / m) * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
    }

    # Combine all probabilities
    all_probs <- c(cluster_probs, aux_probs)

    # Normalize probabilities
    all_probs[is.na(all_probs) | is.infinite(all_probs)] <- 0
    if (sum(all_probs) > 0) {
      all_probs <- all_probs / sum(all_probs)
    } else {
      # Fallback to uniform if all probabilities are 0
      all_probs <- rep(1 / length(all_probs), length(all_probs))
    }

    # Sample new label
    newLabel <- sample.int(length(all_probs), 1, prob = all_probs)

    # Handle the assignment
    if (newLabel <= numLabels) {
      # Assigned to existing cluster
      clusterLabels[i] <- newLabel
      pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    } else {
      # Assigned to auxiliary component - create new cluster
      aux_idx <- newLabel - numLabels

      # Add new cluster
      numLabels <- numLabels + 1
      clusterLabels[i] <- numLabels
      pointsPerCluster <- c(pointsPerCluster, 1)

      # Expand cluster parameters
      for (k in seq_along(clusterParams)) {
        param_dim <- dim(clusterParams[[k]])
        if (length(param_dim) == 3 && param_dim[3] > 0) {
          new_param <- array(NA, dim = c(param_dim[1], param_dim[2], numLabels))
          new_param[, , 1:(numLabels-1)] <- clusterParams[[k]]
          new_param[, , numLabels] <- aux[[k]][, , aux_idx]
          clusterParams[[k]] <- new_param
        } else {
          # Handle case where there were no clusters
          new_param <- array(aux[[k]][, , aux_idx],
                             dim = c(dim(aux[[k]])[1], dim(aux[[k]])[2], 1))
          clusterParams[[k]] <- new_param
        }
      }
    }
  }

  # Final validation - ensure consistency
  if (sum(pointsPerCluster) != n) {
    warning(paste("Inconsistent point counts detected. Expected:", n, "Got:", sum(pointsPerCluster)))
    # Recalculate from cluster labels
    pointsPerCluster <- as.numeric(table(factor(clusterLabels, levels = seq_len(numLabels))))
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
  # Use C++ implementation if enabled and available
  if (using_cpp_hierarchical_samplers() && all(sapply(dpObj$indDP, function(x) inherits(x, "beta")))) {
    return(ClusterComponentUpdate.hierarchical.cpp(dpObj))
  }

  # Original R implementation
  for(i in seq_along(dpObj$indDP)){
    dpObj$indDP[[i]] <- ClusterComponentUpdate(dpObj$indDP[[i]])
    dpObj$indDP[[i]] <- DuplicateClusterRemove(dpObj$indDP[[i]])
  }
  return(dpObj)
}

