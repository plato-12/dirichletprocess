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

  # Check if C++ implementation is available and enabled
  if (inherits(dpObj, "beta") && using_cpp_samplers()) {
    cpp_result <- nonconjugate_beta_cluster_component_update_cpp(dpObj)

    if (!is.null(cpp_result) && !isTRUE(cpp_result$stub_result)) {
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

  # Create working copies
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

    # If cluster is now empty, we'll handle it after assignment
    empty_cluster <- (pointsPerCluster[currentLabel] == 0)

    # Calculate probabilities for existing clusters
    cluster_probs <- numeric(numLabels)

    for (j in 1:numLabels) {
      if (pointsPerCluster[j] > 0 || j == currentLabel) {
        # Extract parameters for cluster j
        if (inherits(dpObj, "beta")) {
          # Special handling for beta distribution parameters
          # Check if parameters are already arrays with correct dimensions
          if (is.array(clusterParams$mu) && length(dim(clusterParams$mu)) == 3) {
            theta_j <- list(
              mu = array(clusterParams$mu[,,j, drop = FALSE], dim = c(1,1,1)),
              nu = array(clusterParams$nu[,,j, drop = FALSE], dim = c(1,1,1))
            )
          } else if (is.array(clusterParams$mu)) {
            # Parameters might be 1D or 2D arrays
            if (length(clusterParams$mu) >= j) {
              theta_j <- list(
                mu = array(clusterParams$mu[j], dim = c(1,1,1)),
                nu = array(clusterParams$nu[j], dim = c(1,1,1))
              )
            } else {
              theta_j <- list(
                mu = array(0.5, dim = c(1,1,1)),  # Default value
                nu = array(1, dim = c(1,1,1))     # Default value
              )
            }
          } else {
            # Parameters are likely vectors or single values
            theta_j <- list(
              mu = array(clusterParams$mu[[j]], dim = c(1,1,1)),
              nu = array(clusterParams$nu[[j]], dim = c(1,1,1))
            )
          }
        } else {
          # Generic parameter extraction for other distributions
          theta_j <- lapply(clusterParams, function(param) {
            if (is.array(param) && length(dim(param)) == 3) {
              array(param[,,j, drop = FALSE],
                    dim = c(dim(param)[1], dim(param)[2], 1))
            } else if (is.array(param) || is.vector(param)) {
              param[j]
            } else if (is.list(param)) {
              param[[j]]
            } else {
              param
            }
          })
        }

        # Calculate likelihood
        lik <- Likelihood(mdObj, y[i, , drop = FALSE], theta_j)

        # Use actual count (0 if empty) for probability calculation
        cluster_probs[j] <- pointsPerCluster[j] * lik
      }
    }

    # Calculate probabilities for auxiliary components
    aux_probs <- numeric(m)
    for (j in 1:m) {
      # Ensure aux parameters have correct structure
      if (inherits(dpObj, "beta")) {
        # Check structure of aux parameters
        if (is.list(dpObj$aux[[j]])) {
          if (is.array(dpObj$aux[[j]]$mu) && length(dim(dpObj$aux[[j]]$mu)) == 3) {
            tempAux <- list(
              mu = array(dpObj$aux[[j]]$mu[,,1, drop = FALSE], dim = c(1,1,1)),
              nu = array(dpObj$aux[[j]]$nu[,,1, drop = FALSE], dim = c(1,1,1))
            )
          } else {
            tempAux <- list(
              mu = array(as.numeric(dpObj$aux[[j]]$mu), dim = c(1,1,1)),
              nu = array(as.numeric(dpObj$aux[[j]]$nu), dim = c(1,1,1))
            )
          }
        } else {
          tempAux <- dpObj$aux[[j]]
        }
      } else {
        tempAux <- lapply(dpObj$aux[[j]], function(param) {
          if (is.array(param) && length(dim(param)) == 3) {
            array(param[,,1, drop = FALSE],
                  dim = c(dim(param)[1], dim(param)[2], 1))
          } else {
            param
          }
        })
      }

      aux_probs[j] <- (alpha / m) * Likelihood(mdObj, y[i, , drop = FALSE], tempAux)
    }

    # Combine all probabilities
    all_probs <- c(cluster_probs, aux_probs)

    # Handle numerical issues
    all_probs[is.na(all_probs) | is.infinite(all_probs) | all_probs < 0] <- 0

    # Normalize probabilities
    prob_sum <- sum(all_probs)
    if (prob_sum > 0) {
      all_probs <- all_probs / prob_sum
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

      # Handle empty cluster removal if necessary
      if (empty_cluster && currentLabel != newLabel) {
        # Remove the empty cluster
        keep_idx <- setdiff(1:numLabels, currentLabel)

        # Create label mapping
        label_map <- integer(numLabels)
        label_map[keep_idx] <- seq_along(keep_idx)

        # Remap all cluster labels
        for (j in seq_len(n)) {
          clusterLabels[j] <- label_map[clusterLabels[j]]
        }

        # Update cluster count
        numLabels <- length(keep_idx)

        # Update points per cluster
        pointsPerCluster <- pointsPerCluster[keep_idx]

        # Update cluster parameters
        if (inherits(dpObj, "beta")) {
          # Special handling for beta parameters
          if (is.array(clusterParams$mu) && length(dim(clusterParams$mu)) == 3) {
            clusterParams$mu <- array(clusterParams$mu[,,keep_idx, drop = FALSE],
                                      dim = c(1, 1, numLabels))
            clusterParams$nu <- array(clusterParams$nu[,,keep_idx, drop = FALSE],
                                      dim = c(1, 1, numLabels))
          } else {
            # Handle other formats
            new_mu <- numeric(numLabels)
            new_nu <- numeric(numLabels)
            for (k in seq_along(keep_idx)) {
              if (is.list(clusterParams$mu)) {
                new_mu[k] <- clusterParams$mu[[keep_idx[k]]]
                new_nu[k] <- clusterParams$nu[[keep_idx[k]]]
              } else {
                new_mu[k] <- clusterParams$mu[keep_idx[k]]
                new_nu[k] <- clusterParams$nu[keep_idx[k]]
              }
            }
            clusterParams$mu <- array(new_mu, dim = c(1, 1, numLabels))
            clusterParams$nu <- array(new_nu, dim = c(1, 1, numLabels))
          }
        } else {
          # Generic parameter update for other distributions
          for (k in seq_along(clusterParams)) {
            param <- clusterParams[[k]]
            if (is.array(param) && length(dim(param)) == 3) {
              clusterParams[[k]] <- array(param[,,keep_idx, drop = FALSE],
                                          dim = c(dim(param)[1], dim(param)[2], numLabels))
            } else if (is.list(param)) {
              clusterParams[[k]] <- param[keep_idx]
            } else if (is.vector(param)) {
              clusterParams[[k]] <- param[keep_idx]
            }
          }
        }
      }

    } else {
      # Assigned to auxiliary component - create new cluster
      aux_idx <- newLabel - numLabels

      if (empty_cluster) {
        # Reuse the empty cluster slot
        clusterLabels[i] <- currentLabel
        pointsPerCluster[currentLabel] <- 1

        # Copy auxiliary parameters to the empty slot
        if (inherits(dpObj, "beta")) {
          if (is.array(clusterParams$mu) && length(dim(clusterParams$mu)) == 3) {
            if (is.array(dpObj$aux[[aux_idx]]$mu) && length(dim(dpObj$aux[[aux_idx]]$mu)) == 3) {
              clusterParams$mu[,,currentLabel] <- dpObj$aux[[aux_idx]]$mu[,,1]
              clusterParams$nu[,,currentLabel] <- dpObj$aux[[aux_idx]]$nu[,,1]
            } else {
              clusterParams$mu[,,currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$mu)
              clusterParams$nu[,,currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$nu)
            }
          } else {
            # Need to maintain structure
            if (is.list(clusterParams$mu)) {
              clusterParams$mu[[currentLabel]] <- as.numeric(dpObj$aux[[aux_idx]]$mu)
              clusterParams$nu[[currentLabel]] <- as.numeric(dpObj$aux[[aux_idx]]$nu)
            } else {
              clusterParams$mu[currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$mu)
              clusterParams$nu[currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$nu)
            }
          }
        } else {
          # Generic parameter copy
          for (k in seq_along(clusterParams)) {
            if (is.array(clusterParams[[k]]) && length(dim(clusterParams[[k]])) == 3) {
              clusterParams[[k]][,,currentLabel] <- dpObj$aux[[aux_idx]][[k]]
            } else {
              clusterParams[[k]][currentLabel] <- dpObj$aux[[aux_idx]][[k]]
            }
          }
        }

      } else {
        # Add new cluster
        numLabels <- numLabels + 1
        clusterLabels[i] <- numLabels
        pointsPerCluster <- c(pointsPerCluster, 1)

        # Expand cluster parameters
        if (inherits(dpObj, "beta")) {
          # Special handling for beta parameters
          if (is.array(clusterParams$mu) && length(dim(clusterParams$mu)) == 3) {
            # Expand 3D arrays
            new_mu <- array(NA, dim = c(1, 1, numLabels))
            new_nu <- array(NA, dim = c(1, 1, numLabels))

            if (numLabels > 1) {
              new_mu[,,1:(numLabels-1)] <- clusterParams$mu
              new_nu[,,1:(numLabels-1)] <- clusterParams$nu
            }

            # Add auxiliary parameter
            if (is.array(dpObj$aux[[aux_idx]]$mu) && length(dim(dpObj$aux[[aux_idx]]$mu)) == 3) {
              new_mu[,,numLabels] <- dpObj$aux[[aux_idx]]$mu[,,1]
              new_nu[,,numLabels] <- dpObj$aux[[aux_idx]]$nu[,,1]
            } else {
              new_mu[,,numLabels] <- as.numeric(dpObj$aux[[aux_idx]]$mu)
              new_nu[,,numLabels] <- as.numeric(dpObj$aux[[aux_idx]]$nu)
            }

            clusterParams$mu <- new_mu
            clusterParams$nu <- new_nu
          } else {
            # Convert to proper 3D structure
            old_mu <- if (is.list(clusterParams$mu)) unlist(clusterParams$mu) else clusterParams$mu
            old_nu <- if (is.list(clusterParams$nu)) unlist(clusterParams$nu) else clusterParams$nu

            new_mu <- array(NA, dim = c(1, 1, numLabels))
            new_nu <- array(NA, dim = c(1, 1, numLabels))

            if (numLabels > 1) {
              for (k in 1:(numLabels-1)) {
                new_mu[,,k] <- old_mu[k]
                new_nu[,,k] <- old_nu[k]
              }
            }

            new_mu[,,numLabels] <- as.numeric(dpObj$aux[[aux_idx]]$mu)
            new_nu[,,numLabels] <- as.numeric(dpObj$aux[[aux_idx]]$nu)

            clusterParams$mu <- new_mu
            clusterParams$nu <- new_nu
          }
        } else {
          # Generic parameter expansion
          for (k in seq_along(clusterParams)) {
            param <- clusterParams[[k]]
            aux_param <- dpObj$aux[[aux_idx]][[k]]

            if (is.array(param) && length(dim(param)) == 3) {
              new_param <- array(NA, dim = c(dim(param)[1], dim(param)[2], numLabels))
              if (numLabels > 1) {
                new_param[,,1:(numLabels-1)] <- param
              }
              new_param[,,numLabels] <- aux_param
              clusterParams[[k]] <- new_param
            } else if (is.list(param)) {
              clusterParams[[k]] <- c(param, list(aux_param))
            } else if (is.vector(param)) {
              clusterParams[[k]] <- c(param, aux_param)
            }
          }
        }
      }
    }
  }

  # Final validation - ensure consistency
  if (sum(pointsPerCluster) != n) {
    warning(paste("Inconsistent point counts detected. Expected:", n,
                  "Got:", sum(pointsPerCluster),
                  "- Recalculating from cluster labels"))

    # Recalculate from cluster labels
    pointsPerCluster <- as.numeric(table(factor(clusterLabels, levels = seq_len(numLabels))))
  }

  # Additional validation checks
  if (any(clusterLabels <= 0)) {
    stop("Non-positive cluster labels detected")
  }

  if (any(clusterLabels > numLabels)) {
    stop(paste("Cluster labels exceed number of clusters.",
               "Max label:", max(clusterLabels),
               "Number of clusters:", numLabels))
  }

  if (length(pointsPerCluster) != numLabels) {
    stop(paste("Length mismatch: pointsPerCluster has", length(pointsPerCluster),
               "elements but numberClusters is", numLabels))
  }

  # Ensure parameters have correct structure for beta
  if (inherits(dpObj, "beta")) {
    if (!is.array(clusterParams$mu) || length(dim(clusterParams$mu)) != 3) {
      # Convert to proper 3D array structure
      mu_vals <- if (is.list(clusterParams$mu)) unlist(clusterParams$mu) else clusterParams$mu
      nu_vals <- if (is.list(clusterParams$nu)) unlist(clusterParams$nu) else clusterParams$nu

      clusterParams$mu <- array(mu_vals, dim = c(1, 1, numLabels))
      clusterParams$nu <- array(nu_vals, dim = c(1, 1, numLabels))
    }
  }

  # Update dpObj with final state
  dpObj$pointsPerCluster <- pointsPerCluster
  dpObj$clusterLabels <- clusterLabels
  dpObj$clusterParameters <- clusterParams
  dpObj$numberClusters <- numLabels

  # Regenerate auxiliary parameters for next iteration
  dpObj$aux <- vector("list", dpObj$m)
  for (j in seq_len(dpObj$m)) {
    dpObj$aux[[j]] <- PriorDraw(dpObj$mixingDistribution, 1)
  }

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

