#' Predict the cluster labels of some new data.
#'
#' Given a fitted Dirichlet process object and some new data use this function to predict what clusters the new data belong to and associated cluster parameters.
#'
#' @param dpobj Fitted Dirichlet Process
#' @param newData New data to have cluster labels predicted.
#' @return A list of the predicted cluster labels of some new unseen data.
#'
#' @examples
#' y <- rnorm(10)
#' dp <- DirichletProcessGaussian(y)
#' dp <- Fit(dp, 5)
#' newY <- rnorm(10, 1)
#' pred <- ClusterLabelPredict(dp, newY)
#'
#' @export
ClusterLabelPredict <- function(dpobj, newData){
  UseMethod("ClusterLabelPredict", dpobj)
}

#' @export
ClusterLabelPredict.conjugate <- function(dpobj, newData) {

  if (!is.matrix(newData))
    newData <- matrix(newData, ncol = 1)

  alpha <- dpobj$alpha
  clusterParams <- dpobj$clusterParameters
  numLabels <- dpobj$numberClusters
  mdobj <- dpobj$mixingDistribution
  pointsPerCluster <- dpobj$pointsPerCluster
  Predictive_newData <- Predictive(mdobj, newData)
  componentIndexes <- numeric(nrow(newData))

  # For mvnormal with pre-allocated arrays, check capacity
  if (inherits(dpobj, "mvnormal") && is.list(clusterParams)) {
    current_capacity <- dim(clusterParams[[1]])[3]
    if (current_capacity < numLabels + nrow(newData)) {
      # Expand arrays preemptively
      new_capacity <- numLabels + nrow(newData) + 20
      for (j in seq_along(clusterParams)) {
        param_dims <- dim(clusterParams[[j]])
        if (length(param_dims) == 3) {
          new_array <- array(NA_real_, dim = c(param_dims[1], param_dims[2], new_capacity))
          new_array[, , 1:current_capacity] <- clusterParams[[j]]

          # Fill remaining slots with prior draws
          if (current_capacity < new_capacity) {
            extra_params <- PriorDraw(mdobj, new_capacity - current_capacity)
            if (j == 1) {
              new_array[, , (current_capacity+1):new_capacity] <- extra_params$mu
            } else {
              new_array[, , (current_capacity+1):new_capacity] <- extra_params$sig
            }
          }

          clusterParams[[j]] <- new_array
        }
      }
    }
  }

  for (i in seq_len(nrow(newData))) {
    dataVal <- newData[i, , drop = FALSE]
    weights <- numeric(numLabels + 1)
    weights[1:numLabels] <- pointsPerCluster * Likelihood(mdobj, dataVal, clusterParams)
    weights[numLabels + 1] <- alpha * Predictive_newData[i]

    ind <- numLabels + 1
    component <- sample.int(ind, 1, prob = weights)

    if (component <= numLabels) {
      componentIndexes[i] <- component
      pointsPerCluster[component] <- pointsPerCluster[component] + 1
    } else {
      componentIndexes[i] <- component
      numLabels <- numLabels + 1
      pointsPerCluster <- c(pointsPerCluster, 1)
      post_draw <- PosteriorDraw(mdobj, newData[i, ,drop=FALSE])

      if (inherits(dpobj, "mvnormal") && is.list(clusterParams)) {
        # For mvnormal with pre-allocated arrays
        current_capacity <- dim(clusterParams[[1]])[3]
        if (numLabels > current_capacity) {
          # This shouldn't happen if we pre-expanded correctly
          stop("Insufficient capacity in pre-allocated arrays")
        }

        # Use the pre-allocated slot
        for (j in seq_along(clusterParams)) {
          param_dims <- dim(clusterParams[[j]])
          if (length(param_dims) == 3) {
            clusterParams[[j]][, , numLabels] <- post_draw[[j]][, , 1]
          }
        }
      } else {
        # Original expansion logic for other distributions
        for (j in seq_along(clusterParams)) {
          clusterParams[[j]] <- array(c(clusterParams[[j]], post_draw[[j]]),
                                      dim = c(dim(post_draw[[j]])[1:2],
                                              dim(clusterParams[[j]])[3] + 1))
        }
      }
    }
  }

  outList <- list(componentIndexes = componentIndexes,
                  pointsPerCluster = pointsPerCluster,
                  clusterParams = clusterParams,
                  numLabels = numLabels)
  return(outList)
}

#' @export
ClusterLabelPredict.nonconjugate <- function(dpobj, newData) {

  if (!is.matrix(newData))
    newData <- matrix(newData, ncol = 1)

  alpha <- dpobj$alpha

  # clusterLabels <- dpobj$clusterLabels
  clusterParams <- dpobj$clusterParameters
  numLabels <- dpobj$numberClusters
  mdobj <- dpobj$mixingDistribution
  m <- dpobj$m

  pointsPerCluster <- dpobj$pointsPerCluster

  componentIndexes <- numeric(length(newData))

  for (i in seq_along(newData)) {

    aux <- PriorDraw(mdobj, m)

    dataVal <- newData[i, , drop = FALSE]
    weights <- numeric(numLabels + 1)

    weights[1:numLabels] <- pointsPerCluster * Likelihood(mdobj, dataVal, clusterParams)
    weights[(numLabels + 1):(numLabels + m)] <- (alpha/m) * Likelihood(mdobj,
      dataVal, aux)

    if (all(weights == 0)) {
      weights[1:(numLabels + m)] <- 1
    }
    if (anyNA(weights)) {
      weights[is.na(weights)] <- 0
    }
    if (any(is.nan(weights))) {
      weights[is.nan(weights)] <- 0
    }

    ind <- numLabels + m
    component <- sample.int(ind, 1, prob = weights)

    if (component <= numLabels) {
      componentIndexes[i] <- component
      pointsPerCluster[component] <- pointsPerCluster[component] + 1
    } else {
      componentIndexes[i] <- numLabels + 1
      pointsPerCluster <- c(pointsPerCluster, 1)
      # clusterParams = rbind(clusterParams, aux[component-numLabels,])

      for (j in seq_along(clusterParams)) {
        clusterParams[[j]] <- array(c(clusterParams[[j]], aux[[j]][, , component - numLabels]), dim = c(dim(clusterParams[[j]])[1:2], dim(clusterParams[[j]])[3] +
          1))
      }

      numLabels <- numLabels + 1
    }
  }
  outList <- list(componentIndexes = componentIndexes,
                  pointsPerCluster = pointsPerCluster,
                  clusterParams = clusterParams,
                  numLabels = numLabels)
  return(outList)
}

