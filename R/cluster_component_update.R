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

    # Generate auxiliary parameters based on whether current cluster is now empty
    if (pointsPerCluster[currentLabel] == 0) {
      # Current cluster is empty, generate m-1 auxiliary parameters
      # and include the current cluster's parameters as one auxiliary
      aux <- PriorDraw(mdObj, m - 1)

      # Combine current cluster params with auxiliary params
      for (k in seq_along(clusterParams)) {
        param_dim <- dim(clusterParams[[k]])
        combined_aux <- array(NA, dim = c(param_dim[1], param_dim[2], m))

        # First auxiliary is the current (empty) cluster's parameters
        combined_aux[, , 1] <- clusterParams[[k]][, , currentLabel]

        # Rest are the newly drawn auxiliary parameters
        if (m > 1) {
          combined_aux[, , 2:m] <- aux[[k]]
        }

        aux[[k]] <- combined_aux
      }
    } else {
      # Current cluster is not empty, generate m auxiliary parameters
      aux <- PriorDraw(mdObj, m)
    }

    # Calculate probabilities for each possible assignment
    probs <- numeric(numLabels + m)

    # Existing clusters
    for (j in seq_len(numLabels)) {
      if (pointsPerCluster[j] > 0) {
        # Non-empty cluster
        tempThetaJ <- vector("list", length(clusterParams))
        for (k in seq_along(clusterParams)) {
          tempThetaJ[[k]] <- array(clusterParams[[k]][, , j],
                                   dim = c(dim(clusterParams[[k]])[1],
                                           dim(clusterParams[[k]])[2], 1))
        }
        probs[j] <- pointsPerCluster[j] * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
      } else {
        # Empty cluster (only possible for currentLabel)
        probs[j] <- 0
      }
    }

    # Auxiliary components
    for (j in seq_len(m)) {
      tempThetaJ <- vector("list", length(clusterParams))
      for (k in seq_along(clusterParams)) {
        tempThetaJ[[k]] <- array(aux[[k]][, , j],
                                 dim = c(dim(aux[[k]])[1],
                                         dim(aux[[k]])[2], 1))
      }
      probs[numLabels + j] <- (alpha / m) * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
    }

    # Normalize probabilities
    probs[is.na(probs)] <- 0
    probs[probs < 0] <- 0

    if (sum(probs) > 0) {
      probs <- probs / sum(probs)
    } else {
      # Fallback to uniform if all probabilities are 0
      probs <- rep(1 / length(probs), length(probs))
    }

    # Sample new label
    newLabel <- sample.int(length(probs), 1, prob = probs)

    # Handle the assignment
    if (newLabel <= numLabels) {
      # Assigned to existing cluster
      clusterLabels[i] <- newLabel
      pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1

      # Check if old cluster is now empty and needs removal
      if (pointsPerCluster[currentLabel] == 0 && currentLabel != newLabel) {
        # Remove empty cluster and shift indices
        numLabels <- numLabels - 1

        # Remove from pointsPerCluster
        if (currentLabel <= length(pointsPerCluster)) {
          pointsPerCluster <- pointsPerCluster[-currentLabel]
        }

        # Remove from cluster parameters
        for (k in seq_along(clusterParams)) {
          param_dim <- dim(clusterParams[[k]])
          if (param_dim[3] > 1) {
            clusterParams[[k]] <- clusterParams[[k]][, , -currentLabel, drop = FALSE]
          }
        }

        # Shift down all labels greater than currentLabel
        clusterLabels[clusterLabels > currentLabel] <- clusterLabels[clusterLabels > currentLabel] - 1

        # Adjust newLabel if it was shifted
        if (newLabel > currentLabel) {
          clusterLabels[i] <- newLabel - 1
        }
      }
    } else {
      # Assigned to auxiliary component - create new cluster
      aux_idx <- newLabel - numLabels

      if (pointsPerCluster[currentLabel] == 0) {
        # Reuse the empty cluster slot
        clusterLabels[i] <- currentLabel
        pointsPerCluster[currentLabel] <- 1

        # Update parameters with auxiliary values
        for (k in seq_along(clusterParams)) {
          clusterParams[[k]][, , currentLabel] <- aux[[k]][, , aux_idx]
        }
      } else {
        # Create genuinely new cluster
        numLabels <- numLabels + 1
        clusterLabels[i] <- numLabels
        pointsPerCluster <- c(pointsPerCluster, 1)

        # Expand cluster parameters
        for (k in seq_along(clusterParams)) {
          param_dim <- dim(clusterParams[[k]])
          new_param <- array(NA, dim = c(param_dim[1], param_dim[2], numLabels))
          new_param[, , 1:(numLabels-1)] <- clusterParams[[k]]
          new_param[, , numLabels] <- aux[[k]][, , aux_idx]
          clusterParams[[k]] <- new_param
        }
      }
    }
  }

  # Final validation - ensure consistency
  if (sum(pointsPerCluster) != n) {
    warning("Inconsistent point counts detected, recalculating...")
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

