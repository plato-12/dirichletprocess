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

    # Compute probabilities for existing clusters - element by element to avoid array issues
    cluster_probs <- numeric(numLabels)

    for (j in 1:numLabels) {
      if (pointsPerCluster[j] > 0) {
        # Create a single-cluster parameter list for likelihood computation
        single_cluster_params <- list(
          array(clusterParams[[1]][, , j], dim = c(1, 1, 1)),
          array(clusterParams[[2]][, , j], dim = c(1, 1, 1))
        )

        # Compute likelihood for this specific cluster
        likelihood_val <- Likelihood(mdObj, y[i, , drop = FALSE], single_cluster_params)
        cluster_probs[j] <- pointsPerCluster[j] * as.numeric(likelihood_val[1])
      } else {
        cluster_probs[j] <- 0
      }
    }

    # Probability for new cluster
    new_cluster_prob <- alpha * predictiveArray[i]

    # Combine all probabilities
    probs <- c(cluster_probs, new_cluster_prob)

    # Handle edge cases
    probs[is.na(probs) | is.infinite(probs)] <- 0

    if (all(probs == 0)) {
      probs <- rep_len(1, length(probs))
    }

    # Sample new cluster assignment
    newLabel <- sample.int(numLabels + 1, 1, prob = probs)

    # Restore the point count before calling ClusterLabelChange
    pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] + 1
    dpObj$pointsPerCluster <- pointsPerCluster

    # Apply the cluster change
    dpObj <- ClusterLabelChange(dpObj, i, newLabel, currentLabel)

    # Update local variables from the modified dpObj
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

    if (pointsPerCluster[currentLabel] == 0) {

      priorDraws <- PriorDraw(mdObj, m - 1)

      for (j in seq_along(priorDraws)) {
        aux[[j]] <- array(c(clusterParams[[j]][, , currentLabel], priorDraws[[j]]),
                          dim = c(dim(priorDraws[[j]])[1:2], m))
      }
    } else {
      aux <- PriorDraw(mdObj, m)
    }

    probs <- c(
      pointsPerCluster * Likelihood(mdObj, y[i, , drop = FALSE],clusterParams),
      (alpha/m) * Likelihood(mdObj, y[i, , drop = FALSE], aux))

    if (any(is.nan(probs))) {
      probs[is.nan(probs)] <- 0
    }


    probs[is.na(probs)] <- 0


    if (any(is.infinite(probs))) {
      probs[is.infinite(probs)] <- 1
      probs[-is.infinite(probs)] <- 0
    }

    if (all(probs == 0)) {
      probs <- rep_len(1, length(probs))
    }
    newLabel <- sample.int(numLabels + m, 1, prob = probs)

    dpObj$pointsPerCluster <- pointsPerCluster

    dpObj <- ClusterLabelChange(dpObj, i, newLabel, currentLabel, aux)

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

#' @export
#' @rdname ClusterComponentUpdate
ClusterComponentUpdate.hierarchical <- function(dpObj){

  for(i in seq_along(dpObj$indDP)){
    dpObj$indDP[[i]] <- ClusterComponentUpdate(dpObj$indDP[[i]])
    dpObj$indDP[[i]] <- DuplicateClusterRemove(dpObj$indDP[[i]])
  }
  return(dpObj)
}
