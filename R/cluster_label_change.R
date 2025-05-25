#' Change cluster labels in a Dirichlet Process object
#'
#' Internal function to handle cluster label changes, including creation of new clusters
#' and removal of empty clusters.
#'
#' @param dpObj Dirichlet process object
#' @param i Index of the data point to reassign
#' @param newLabel New cluster label for the data point
#' @param currentLabel Current cluster label of the data point
#' @param aux Auxiliary parameters for non-conjugate case
#' @return Updated Dirichlet process object
#' @export
ClusterLabelChange <- function(dpObj, i, newLabel, currentLabel, aux=0) {
  UseMethod("ClusterLabelChange", dpObj)
}

#' @export
ClusterLabelChange.conjugate <- function(dpObj, i, newLabel, currentLabel, aux=0) {

  x <- dpObj$data[i, , drop = FALSE]
  pointsPerCluster <- dpObj$pointsPerCluster
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  mdObj <- dpObj$mixingDistribution

  # Caller of ClusterLabelChange is responsible for decrementing pointsPerCluster[currentLabel]
  # This function handles the assignment to newLabel and potential cleanup/creation.

  if (newLabel <= numLabels) { # Assigning to an existing cluster slot
    pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    clusterLabels[i] <- newLabel

    # If the original cluster (currentLabel) is now empty (caller would have decremented it,
    # and if it was the only point, it's now 0 AFTER that decrement by caller)
    # This check should be pointsPerCluster[currentLabel] == 0
    # (assuming the caller already removed the point from currentLabel for counting purposes)
    # However, standard algorithm 8 implies this function does the full management.
    # Let's assume currentLabel's count includes point i *before* this function's logic fully applied the change.
    # The critical part is what the caller (ClusterComponentUpdate) does *before* calling this.
    # If caller does: dpObj$pointsPerCluster[ci] <- dpObj$pointsPerCluster[ci] - 1
    # then currentLabel's count is already reduced.

    if (pointsPerCluster[currentLabel] == 0 && currentLabel != newLabel) { # Old cluster (not the one we moved to) is now empty
      numLabels <- numLabels - 1
      pointsPerCluster <- pointsPerCluster[-currentLabel]
      clusterParams <- lapply(clusterParams, function(param_array) param_array[, , -currentLabel, drop = FALSE])

      # Adjust labels for clusters that were after the removed one
      # And also adjust newLabel if it was affected by the shift
      original_newLabel_val <- newLabel # Store before potential shift

      inds_labels_to_decrement <- clusterLabels > currentLabel
      clusterLabels[inds_labels_to_decrement] <- clusterLabels[inds_labels_to_decrement] - 1

      # If newLabel itself was shifted down because it was > currentLabel
      if (original_newLabel_val > currentLabel) {
        clusterLabels[i] <- original_newLabel_val - 1 # Correctly point to the shifted newLabel
      }
    }
  } else { # Assigning to a new cluster (newLabel refers to an auxiliary parameter index)
    # newLabel from caller is an index into a combined list of existing + m auxiliary parameters.
    # If it's > numLabels, it means one of the m auxiliary parameters was chosen.
    # The actual new cluster index will be numLabels + 1.

    # Case 1: The point's original cluster (currentLabel) becomes empty
    if (pointsPerCluster[currentLabel] == 0) {
      # Re-purpose the slot of the now-empty currentLabel for the new cluster parameters
      post_draw <- PosteriorDraw(mdObj, x) # Parameters for the new cluster based on point x
      for (k in seq_along(clusterParams)) {
        clusterParams[[k]][, , currentLabel] <- post_draw[[k]]
      }
      pointsPerCluster[currentLabel] <- 1 # This slot now has point i
      clusterLabels[i] <- currentLabel   # Point i is assigned to this re-purposed cluster index
      # numLabels does not change because we reused a slot.
    } else {
      # Case 2: The point's original cluster is not empty, so we truly add a new cluster
      clusterLabels[i] <- numLabels + 1 # Assign to the next available cluster index

      pointsPerCluster <- c(pointsPerCluster, 1) # Add count for the new cluster

      post_draw <- PosteriorDraw(mdObj, x) # Parameters for the new cluster
      for (j in seq_along(clusterParams)) {
        dim_existing <- dim(clusterParams[[j]])
        new_param_array <- array(NA, dim = c(dim_existing[1], dim_existing[2], numLabels + 1)) # new size
        if(numLabels > 0 && dim_existing[3] > 0) { # Copy if there were existing params
          new_param_array[,,1:numLabels] <- clusterParams[[j]]
        }
        new_param_array[,,numLabels+1] <- post_draw[[j]]
        clusterParams[[j]] <- new_param_array
      }
      numLabels <- numLabels + 1 # Increment total number of clusters
    }
  }

  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels
  return(dpObj)
}

#' @export
ClusterLabelChange.nonconjugate <- function(dpObj, i, newLabel, currentLabel, aux=0) {

  pointsPerCluster <- dpObj$pointsPerCluster
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters # numLabels before potential change

  # As in conjugate, assume caller (ClusterComponentUpdate) has already decremented pointsPerCluster[currentLabel]

  if (newLabel <= numLabels) { # Assigning to an existing, non-empty cluster (or an empty one that's not currentLabel)
    pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    clusterLabels[i] <- newLabel

    if (pointsPerCluster[currentLabel] == 0 && currentLabel != newLabel) { # Old cluster is now empty and it's not the one we moved to
      numLabels <- numLabels - 1
      pointsPerCluster <- pointsPerCluster[-currentLabel]
      clusterParams <- lapply(clusterParams, function(x) x[, , -currentLabel, drop = FALSE])

      original_newLabel_val <- newLabel
      inds <- clusterLabels > currentLabel
      clusterLabels[inds] <- clusterLabels[inds] - 1
      if (original_newLabel_val > currentLabel) {
        clusterLabels[i] <- original_newLabel_val - 1
      }
    }
  } else { # Assigning to what will become a new cluster, using one of the auxiliary parameters
    aux_param_index_in_aux_list <- newLabel - numLabels # Index into the 'aux' list of parameters (1 to m)

    if (pointsPerCluster[currentLabel] == 0) { # Old cluster became empty, reuse its slot for the new cluster
      for (j in seq_along(clusterParams)) {
        clusterParams[[j]][, , currentLabel] <- aux[[j]][, , aux_param_index_in_aux_list, drop = FALSE]
      }
      pointsPerCluster[currentLabel] <- 1
      clusterLabels[i] <- currentLabel
      # numLabels does not change
    } else { # Old cluster not empty, truly adding a new cluster
      clusterLabels[i] <- numLabels + 1
      pointsPerCluster <- c(pointsPerCluster, 1)

      for (j in seq_along(clusterParams)) {
        dim_existing <- dim(clusterParams[[j]])
        new_param_array <- array(NA, dim = c(dim_existing[1], dim_existing[2], numLabels + 1))
        if(numLabels > 0 && dim_existing[3] > 0) {
          new_param_array[,,1:numLabels] <- clusterParams[[j]]
        }
        new_param_array[,,numLabels+1] <- aux[[j]][, , aux_param_index_in_aux_list, drop = FALSE]
        clusterParams[[j]] <- new_param_array
      }
      numLabels <- numLabels + 1
    }
  }

  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels
  return(dpObj)
}

#' @export
ClusterLabelChange.default <- function(dpObj, i, newLabel, currentLabel, aux=0) {
  # All cat() and print() statements specific to "R-CLC" were previously removed.
  # The "R-CCU:" prints are originating from elsewhere (likely test files or ClusterComponentUpdate.R).
  if (inherits(dpObj, "conjugate")) {
    return(ClusterLabelChange.conjugate(dpObj, i, newLabel, currentLabel, aux))
  } else if (inherits(dpObj, "nonconjugate")) {
    return(ClusterLabelChange.nonconjugate(dpObj, i, newLabel, currentLabel, aux))
  } else {
    stop("ClusterLabelChange not implemented for this object type")
  }
}
