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
        mu = array(clusterParams[[1]][, , currentLabel], dim = c(1, dim(clusterParams[[1]])[2], 1)),
        sig = array(clusterParams[[2]][, , currentLabel], dim = c(dim(clusterParams[[2]])[1], dim(clusterParams[[2]])[2], 1))
      )
    } else {
      NULL
    }

    if (!is.null(current_params_for_empty_slot)) {
      # Get m-1 auxiliary parameters
      aux_params <- PriorDraw(mdObj, m - 1)

      # Combine current cluster params with auxiliary parameters
      for (j in seq_along(clusterParams)) {
        param_dim <- dim(clusterParams[[j]])
        if (length(param_dim) == 3) {
          # Create array to hold all m parameters
          new_dim <- param_dim
          new_dim[3] <- m
          aux[[j]] <- array(NA, dim = new_dim)

          # First position gets the current (soon to be empty) cluster's parameters
          if (length(param_dim) == 3) {
            aux[[j]][, , 1] <- clusterParams[[j]][, , currentLabel]
          } else {
            aux[[j]][, , 1] <- current_params_for_empty_slot[[j]][, , 1]
          }

          # Rest get the auxiliary parameters
          for (k in 2:m) {
            aux[[j]][, , k] <- aux_params[[j]][, , k-1]
          }
        }
      }
    } else {
      aux <- PriorDraw(mdObj, m)
    }

    probs <- numeric(numLabels + m)

    for (j in seq_len(numLabels)) {
      if (j == currentLabel && pointsPerCluster[j] == 0) {
        # Use auxiliary parameter for empty cluster
        tempThetaJ <- vector("list", length(clusterParams))
        for (k in seq_along(clusterParams)) {
          tempThetaJ[[k]] <- array(aux[[k]][, , 1], dim = c(dim(aux[[k]])[1], dim(aux[[k]])[2], 1))
        }

        probs[j] <- (1) * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
      } else if (pointsPerCluster[j] > 0) {
        tempThetaJ <- vector("list", length(clusterParams))
        for (k in seq_along(clusterParams)) {
          tempThetaJ[[k]] <- array(clusterParams[[k]][, , j],
                                   dim = c(dim(clusterParams[[k]])[1], dim(clusterParams[[k]])[2], 1))
        }

        probs[j] <- pointsPerCluster[j] * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
      }
    }

    # Auxiliary components
    for (j in seq_len(m)) {
      tempThetaJ <- vector("list", length(clusterParams))
      for (k in seq_along(clusterParams)) {
        if (pointsPerCluster[currentLabel] == 0 && j == 1) {
          # Skip the first aux parameter if current cluster is empty (already used above)
          next
        }
        aux_idx <- if (pointsPerCluster[currentLabel] == 0 && j > 1) j else j
        tempThetaJ[[k]] <- array(aux[[k]][, , aux_idx],
                                 dim = c(dim(aux[[k]])[1], dim(aux[[k]])[2], 1))
      }

      probs[numLabels + j] <- (alpha / m) * Likelihood(mdObj, y[i, , drop = FALSE], tempThetaJ)
    }

    probs[is.na(probs)] <- 0
    probs[probs < 0] <- 0

    if (sum(probs) > 0) {
      probs <- probs / sum(probs)
    } else {
      probs <- rep(1 / length(probs), length(probs))
    }

    newLabel <- sample.int(length(probs), 1, prob = probs)

    # Update the cluster assignment
    if (newLabel <= numLabels) {
      # Assigned to existing cluster
      clusterLabels[i] <- newLabel
      pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    } else {
      # Assigned to auxiliary component - create new cluster
      numLabels <- numLabels + 1
      clusterLabels[i] <- numLabels
      pointsPerCluster[numLabels] <- 1

      aux_component_index <- newLabel - (numLabels - 1)

      # Expand cluster parameters
      for (j in seq_along(clusterParams)) {
        param_dim <- dim(clusterParams[[j]])
        new_dim <- param_dim
        new_dim[3] <- numLabels

        new_param <- array(NA, dim = new_dim)
        if (param_dim[3] > 0) {
          new_param[, , 1:param_dim[3]] <- clusterParams[[j]]
        }
        new_param[, , numLabels] <- aux[[j]][, , aux_component_index]
        clusterParams[[j]] <- new_param
      }
    }

    # Handle empty clusters by removing them
    if (pointsPerCluster[currentLabel] == 0 && currentLabel <= numLabels) {
      # Remove empty cluster
      if (currentLabel == numLabels) {
        # It's the last cluster, just reduce count
        numLabels <- numLabels - 1

        # Shrink parameter arrays
        for (k in seq_along(clusterParams)) {
          param_dim <- dim(clusterParams[[k]])
          if (length(param_dim) == 3 && param_dim[3] > 1) {
            clusterParams[[k]] <- clusterParams[[k]][, , 1:(param_dim[3]-1), drop = FALSE]
          }
        }

        # Remove last element from pointsPerCluster
        pointsPerCluster <- pointsPerCluster[1:numLabels]

      } else {
        # Not the last cluster - need to shift everything down

        # Shift labels down for clusters after the empty one
        clusterLabels[clusterLabels > currentLabel] <- clusterLabels[clusterLabels > currentLabel] - 1

        # Shift parameters
        for (k in seq_along(clusterParams)) {
          param_dim <- dim(clusterParams[[k]])
          if (length(param_dim) == 3 && param_dim[3] >= currentLabel) {
            # Create new array without the empty cluster
            new_dim <- param_dim
            new_dim[3] <- param_dim[3] - 1
            new_param <- array(NA, dim = new_dim)

            # Copy parameters before the empty cluster
            if (currentLabel > 1) {
              new_param[, , 1:(currentLabel-1)] <- clusterParams[[k]][, , 1:(currentLabel-1)]
            }

            # Copy parameters after the empty cluster (shifted down)
            if (currentLabel < param_dim[3]) {
              new_param[, , currentLabel:(new_dim[3])] <- clusterParams[[k]][, , (currentLabel+1):param_dim[3]]
            }

            clusterParams[[k]] <- new_param
          }
        }

        # Shift pointsPerCluster
        new_pointsPerCluster <- numeric(length(pointsPerCluster) - 1)
        if (currentLabel > 1) {
          new_pointsPerCluster[1:(currentLabel-1)] <- pointsPerCluster[1:(currentLabel-1)]
        }
        if (currentLabel < length(pointsPerCluster)) {
          new_pointsPerCluster[currentLabel:(length(new_pointsPerCluster))] <-
            pointsPerCluster[(currentLabel+1):length(pointsPerCluster)]
        }
        pointsPerCluster <- new_pointsPerCluster

        numLabels <- numLabels - 1
      }
    }
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

#' @export
ClusterComponentUpdate.beta <- function(dpObj) {

  # Get auxiliary parameters
  aux <- PriorDraw(dpObj$mixingDistribution, dpObj$m)

  # Store old labels for comparison
  oldLabels <- dpObj$clusterLabels

  for (i in seq_len(dpObj$n)) {

    currentLabel <- dpObj$clusterLabels[i]

    # Remove point from current cluster
    dpObj$pointsPerCluster[currentLabel] <- dpObj$pointsPerCluster[currentLabel] - 1

    if (dpObj$pointsPerCluster[currentLabel] == 0) {
      # Remove empty cluster
      dpObj$clusterParameters[[1]] <- dpObj$clusterParameters[[1]][, , -currentLabel, drop = FALSE]
      dpObj$clusterParameters[[2]] <- dpObj$clusterParameters[[2]][, , -currentLabel, drop = FALSE]

      dpObj$pointsPerCluster <- dpObj$pointsPerCluster[-currentLabel]
      dpObj$numberClusters <- dpObj$numberClusters - 1

      # Adjust labels
      adjustedLabels <- dpObj$clusterLabels
      adjustedLabels[dpObj$clusterLabels > currentLabel] <- adjustedLabels[dpObj$clusterLabels > currentLabel] - 1
      dpObj$clusterLabels <- adjustedLabels
    }

    # Calculate probabilities for existing clusters
    clusterProbs <- numeric(dpObj$numberClusters)

    for (j in seq_len(dpObj$numberClusters)) {
      theta_j <- list(
        mu = dpObj$clusterParameters[[1]][, , j, drop = FALSE],
        nu = dpObj$clusterParameters[[2]][, , j, drop = FALSE]
      )
      clusterProbs[j] <- dpObj$pointsPerCluster[j] *
        Likelihood(dpObj$mixingDistribution, dpObj$data[i], theta_j)
    }

    # Add auxiliary clusters
    auxProbs <- numeric(dpObj$m)
    for (j in seq_len(dpObj$m)) {
      theta_aux <- list(
        mu = aux[[1]][, , j, drop = FALSE],
        nu = aux[[2]][, , j, drop = FALSE]
      )
      auxProbs[j] <- (dpObj$alpha/dpObj$m) *
        Likelihood(dpObj$mixingDistribution, dpObj$data[i], theta_aux)
    }

    # Sample new cluster
    probs <- c(clusterProbs, auxProbs)
    probs <- probs / sum(probs)

    newCluster <- sample.int(length(probs), 1, prob = probs)

    if (newCluster <= dpObj$numberClusters) {
      # Assign to existing cluster
      dpObj$clusterLabels[i] <- newCluster
      dpObj$pointsPerCluster[newCluster] <- dpObj$pointsPerCluster[newCluster] + 1
    } else {
      # Create new cluster
      dpObj$numberClusters <- dpObj$numberClusters + 1
      dpObj$clusterLabels[i] <- dpObj$numberClusters
      dpObj$pointsPerCluster <- c(dpObj$pointsPerCluster, 1)

      # Add parameters from auxiliary
      auxIndex <- newCluster - length(clusterProbs)
      newMu <- aux[[1]][, , auxIndex, drop = FALSE]
      newNu <- aux[[2]][, , auxIndex, drop = FALSE]

      dpObj$clusterParameters[[1]] <- abind(dpObj$clusterParameters[[1]], newMu, along = 3)
      dpObj$clusterParameters[[2]] <- abind(dpObj$clusterParameters[[2]], newNu, along = 3)
    }
  }

  return(dpObj)
}
