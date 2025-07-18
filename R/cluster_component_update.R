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
        # Extract the parameters for cluster j, preserving dimensions and names
        single_cluster_params <- list()
        for (k in seq_along(clusterParams)) {
          param_name <- names(clusterParams)[k]
          if (is.null(param_name) || param_name == "") {
            # If no name, use the index (fallback)
            param_name <- k
          }
          param_dims <- dim(clusterParams[[k]])
          if (length(param_dims) == 3) {
            # For 3D arrays (FULL covariance models), extract the slice for cluster j
            single_cluster_params[[param_name]] <- array(
              clusterParams[[k]][, , j],
              dim = c(param_dims[1], param_dims[2], 1)
            )
          } else if (length(param_dims) == 2) {
            # For 2D arrays (constrained covariance models), extract column j
            single_cluster_params[[param_name]] <- clusterParams[[k]][, j, drop = FALSE]
          } else {
            # Fallback for 1D or scalar structures
            single_cluster_params[[param_name]] <- clusterParams[[k]][j]
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

#' @export
#' @rdname ClusterComponentUpdate
ClusterComponentUpdate.nonconjugate <- function(dpObj) {
  # Use C++ implementation if enabled and available
  if (using_cpp()) {
    if (inherits(dpObj, "beta") && exists("nonconjugate_beta_cluster_component_update_cpp")) {
      return(nonconjugate_beta_cluster_component_update_cpp(dpObj))
    }
    if (inherits(dpObj, "mvnormal2") && exists("nonconjugate_mvnormal2_cluster_component_update_cpp")) {
      return(nonconjugate_mvnormal2_cluster_component_update_cpp(dpObj))
    }
  }

  y <- dpObj$data
  n <- dpObj$n
  alpha <- dpObj$alpha
  m <- dpObj$m

  clusterLabels <- dpObj$clusterLabels
  clusterParams <- dpObj$clusterParameters
  numLabels <- dpObj$numberClusters
  mdObj <- dpObj$mixingDistribution

  pointsPerCluster <- dpObj$pointsPerCluster

  for (i in seq_len(n)) {
    currentLabel <- clusterLabels[i]
    pointsPerCluster[currentLabel] <- pointsPerCluster[currentLabel] - 1

    empty_cluster <- (pointsPerCluster[currentLabel] == 0)

    # Calculate probabilities for existing clusters
    cluster_probs <- numeric(numLabels)
    for (j in seq_len(numLabels)) {
      if (pointsPerCluster[j] > 0 || (empty_cluster && j == currentLabel)) {
        # Extract parameters for cluster j
        single_cluster_params <- list()
        for (k in seq_along(clusterParams)) {
          param_dims <- dim(clusterParams[[k]])
          if (length(param_dims) == 3) {
            single_cluster_params[[k]] <- array(
              clusterParams[[k]][, , j],
              dim = c(param_dims[1], param_dims[2], 1)
            )
          } else {
            single_cluster_params[[k]] <- clusterParams[[k]][j]
          }
        }

        # Calculate likelihood
        lik_val <- Likelihood(mdObj, y[i, , drop = FALSE], single_cluster_params)

        # Handle NA/NaN/Inf values
        if (is.na(lik_val) || is.nan(lik_val) || is.infinite(lik_val)) {
          lik_val <- 0
        } else if (lik_val < 0) {
          lik_val <- 0
        }

        if (empty_cluster && j == currentLabel) {
          cluster_probs[j] <- (alpha / m) * lik_val
        } else {
          cluster_probs[j] <- pointsPerCluster[j] * lik_val
        }
      } else {
        cluster_probs[j] <- 0
      }
    }

    # Calculate probabilities for auxiliary components
    aux_probs <- numeric(m)
    for (j in seq_len(m)) {
      lik_val <- Likelihood(mdObj, y[i, , drop = FALSE], dpObj$aux[[j]])

      # Handle NA/NaN/Inf values
      if (is.na(lik_val) || is.nan(lik_val) || is.infinite(lik_val)) {
        lik_val <- 0
      } else if (lik_val < 0) {
        lik_val <- 0
      }

      aux_probs[j] <- (alpha / m) * lik_val
    }

    # Combine probabilities
    all_probs <- c(cluster_probs, aux_probs)

    # Additional safety check for all zeros or invalid values
    if (all(all_probs <= 0) || all(is.na(all_probs)) || sum(all_probs) == 0) {
      # Fallback to uniform distribution
      all_probs <- rep(1, length(all_probs))
    }

    # Sample new label
    newLabel <- sample.int(numLabels + m, 1, prob = all_probs)

    # Update cluster assignment
    if (newLabel <= numLabels) {
      # Assigned to existing cluster
      clusterLabels[i] <- newLabel
      pointsPerCluster[newLabel] <- pointsPerCluster[newLabel] + 1

      # Clean up empty cluster if needed
      if (empty_cluster && newLabel != currentLabel) {
        # Remove empty cluster
        keep_idx <- seq_len(numLabels)[-currentLabel]

        # Update labels - ensure labels remain positive
        shift_mask <- clusterLabels > currentLabel
        clusterLabels[shift_mask] <- clusterLabels[shift_mask] - 1
        
        # Adjust the current point's label if it was affected by the shift
        if (newLabel > currentLabel) {
          clusterLabels[i] <- newLabel - 1
        }

        # Update number of clusters
        numLabels <- numLabels - 1

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
          # Special handling for beta parameters
          if (is.array(clusterParams$mu) && length(dim(clusterParams$mu)) == 3) {
            if (is.list(dpObj$aux[[aux_idx]]) && all(c("mu", "nu") %in% names(dpObj$aux[[aux_idx]]))) {
              clusterParams$mu[,,currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$mu[,,1])
              clusterParams$nu[,,currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$nu[,,1])
            } else {
              # Generate new parameters if aux structure is wrong
              newParams <- PriorDraw(dpObj$mixingDistribution, 1)
              clusterParams$mu[,,currentLabel] <- as.numeric(newParams$mu)
              clusterParams$nu[,,currentLabel] <- as.numeric(newParams$nu)
            }
          } else {
            # Handle other parameter formats
            if (is.list(dpObj$aux[[aux_idx]])) {
              clusterParams[[1]][currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$mu[,,1])
              clusterParams[[2]][currentLabel] <- as.numeric(dpObj$aux[[aux_idx]]$nu[,,1])
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
            if (is.list(dpObj$aux[[aux_idx]]) && all(c("mu", "nu") %in% names(dpObj$aux[[aux_idx]]))) {
              new_mu[,,numLabels] <- as.numeric(dpObj$aux[[aux_idx]]$mu[,,1])
              new_nu[,,numLabels] <- as.numeric(dpObj$aux[[aux_idx]]$nu[,,1])
            } else {
              # Generate new parameters if aux structure is wrong
              newParams <- PriorDraw(dpObj$mixingDistribution, 1)
              new_mu[,,numLabels] <- as.numeric(newParams$mu)
              new_nu[,,numLabels] <- as.numeric(newParams$nu)
            }

            clusterParams$mu <- new_mu
            clusterParams$nu <- new_nu
          } else {
            # Handle other formats
            clusterParams[[1]] <- c(clusterParams[[1]], as.numeric(dpObj$aux[[aux_idx]]$mu[,,1]))
            clusterParams[[2]] <- c(clusterParams[[2]], as.numeric(dpObj$aux[[aux_idx]]$nu[,,1]))
          }
        } else {
          # Generic expansion for other distributions
          for (k in seq_along(clusterParams)) {
            if (is.array(clusterParams[[k]]) && length(dim(clusterParams[[k]])) == 3) {
              old_dim <- dim(clusterParams[[k]])
              new_array <- array(NA, dim = c(old_dim[1], old_dim[2], numLabels))
              if (numLabels > 1) {
                new_array[,,1:(numLabels-1)] <- clusterParams[[k]]
              }
              new_array[,,numLabels] <- dpObj$aux[[aux_idx]][[k]]
              clusterParams[[k]] <- new_array
            } else {
              clusterParams[[k]] <- c(clusterParams[[k]], dpObj$aux[[aux_idx]][[k]])
            }
          }
        }
      }
    }
  }

  # Additional validation checks - fix cluster labels first
  if (any(clusterLabels <= 0)) {
    # Fix non-positive cluster labels
    min_label <- min(clusterLabels)
    if (min_label <= 0) {
      clusterLabels <- clusterLabels - min_label + 1
      warning("Non-positive cluster labels detected and fixed")
    }
  }

  # Update numLabels based on actual cluster labels
  numLabels <- max(clusterLabels)

  # Final validation - recalculate pointsPerCluster after any label fixes
  if (sum(pointsPerCluster) != n || length(pointsPerCluster) != numLabels) {
    warning(paste("Points per cluster mismatch after update.",
                  "Expected:", n,
                  "Got:", sum(pointsPerCluster),
                  "- Recalculating from cluster labels"))

    # Recalculate from cluster labels
    pointsPerCluster <- as.numeric(table(factor(clusterLabels, levels = seq_len(numLabels))))
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

