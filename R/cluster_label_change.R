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

  if (newLabel <= numLabels) {
    pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    clusterLabels[i] <- newLabel

    if (pointsPerCluster[currentLabel] == 0) {
      ### Removing the Empty Cluster ###
      numLabels <- numLabels - 1
      pointsPerCluster <- pointsPerCluster[-currentLabel]

      # clusterParams <- clusterParams[-currentLabel, ,drop=FALSE]
      clusterParams <- lapply(clusterParams, function(x) x[, , -currentLabel,
                                                           drop = FALSE])

      inds <- clusterLabels > currentLabel
      clusterLabels[inds] <- clusterLabels[inds] - 1
    }
  } else {

    if (pointsPerCluster[currentLabel] == 0) {

      post_draw <- PosteriorDraw(mdObj, x)

      for (i in seq_along(clusterParams)) {
        clusterParams[[i]][, , currentLabel] <- post_draw[[i]]
      }

      pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] + 1
    } else {

      clusterLabels[i] <- newLabel
      numLabels <- numLabels + 1
      pointsPerCluster <- c(pointsPerCluster, 1)

      post_draw <- PosteriorDraw(mdObj, x)

      # clusterParams = rbind(clusterParams, posteriorDraw(mdObj, x))

      for (j in seq_along(clusterParams)) {
        clusterParams[[j]] <- array(c(clusterParams[[j]], post_draw[[j]]),
                                    dim = c(dim(post_draw[[j]])[1:2], dim(clusterParams[[j]])[3] +
                                              1))
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
ClusterLabelChange.nonconjugate <- function(dpObj, i, newLabel, currentLabel, aux=0) {

  pointsPerCluster <- dpObj$pointsPerCluster
  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  # mdObj <- dpObj$mixingDistribution

  if (newLabel <= numLabels) {
    pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1
    clusterLabels[i] <- newLabel

    if (pointsPerCluster[currentLabel] == 0) {
      # print('B') Removing the Empty Cluster ###
      numLabels <- numLabels - 1
      pointsPerCluster <- pointsPerCluster[-currentLabel]
      # clusterParams <- clusterParams[-currentLabel, ,drop=FALSE]
      clusterParams <- lapply(clusterParams, function(x) x[, , -currentLabel,
                                                           drop = FALSE])

      inds <- clusterLabels > currentLabel
      clusterLabels[inds] <- clusterLabels[inds] - 1
    }
  } else {

    if (pointsPerCluster[currentLabel] == 0) {
      # print('C') clusterParams[currentLabel, ] = aux[newLabel-numLabels, ]

      for (j in seq_along(clusterParams)) {
        clusterParams[[j]][, , currentLabel] <- aux[[j]][, , newLabel - numLabels]
      }
      pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] + 1

    } else {
      # print('D')
      clusterLabels[i] <- numLabels + 1
      pointsPerCluster <- c(pointsPerCluster, 1)
      # clusterParams = rbind(clusterParams, aux[newLabel-numLabels, ])

      for (j in seq_along(clusterParams)) {
        clusterParams[[j]] <- array(c(clusterParams[[j]],
                                      aux[[j]][, , newLabel - numLabels]),
                                    dim = c(dim(clusterParams[[j]])[1:2],
                                            dim(clusterParams[[j]])[3] + 1))
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
  # Determine conjugacy and dispatch appropriately
  if (inherits(dpObj, "conjugate")) {
    return(ClusterLabelChange.conjugate(dpObj, i, newLabel, currentLabel, aux))
  } else if (inherits(dpObj, "nonconjugate")) {
    return(ClusterLabelChange.nonconjugate(dpObj, i, newLabel, currentLabel, aux))
  } else {
    stop("ClusterLabelChange not implemented for this object type")
  }
}
