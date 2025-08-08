#' Fit the Dirichlet process object
#'
#' Using Neal's algorithm 4 or 8 depending on conjugacy the sampling procedure
#' for a Dirichlet process is carried out. Lists of both cluster parameters,
#' weights and the sampled concentration values are included in the fitted dpObj.
#' When update_prior is set to TRUE the parameters of the base measure are also updated.
#'
#' @param dpObj Initialised Dirichlet Process object
#' @param its Number of iterations to use
#' @param updatePrior Logical flag, defaults to FALSE. Set whether the parameters
#'        of the base measure are updated.
#' @param progressBar Logical flag indicating whether to display a progress bar.
#' @param ... Additional arguments
#' @return A Dirichlet Process object with the fitted cluster parameters and labels.
#'
#' @references Neal, R. M. (2000). Markov chain sampling methods for Dirichlet
#'             process mixture models. Journal of computational and graphical
#'             statistics, 9(2), 249-265.
#'
#' @export
Fit <- function(dpObj, its, updatePrior = FALSE, progressBar = TRUE, ...) {
  UseMethod("Fit", dpObj)
}

#' @export
Fit.default <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive(), ...) {

  if (progressBar) {
    pb <- txtProgressBar(min = 0, max = its, width = 50, char = "-", style = 3)
  }

  alphaChain <- numeric(its)
  likelihoodChain <- numeric(its)
  weightsChain <- vector("list", length = its)
  clusterParametersChain <- vector("list", length = its)
  priorParametersChain <- vector("list", length = its)
  labelsChain <- vector("list", length = its)

  for (i in seq_len(its)) {

    alphaChain[i] <- dpObj$alpha
    weightsChain[[i]] <- dpObj$pointsPerCluster / dpObj$n
    clusterParametersChain[[i]] <- dpObj$clusterParameters
    priorParametersChain[[i]] <- dpObj$mixingDistribution$priorParameters
    labelsChain[[i]] <- dpObj$clusterLabels

    likelihoodChain[i] <- sum(log(LikelihoodDP(dpObj)))

    dpObj <- ClusterComponentUpdate(dpObj)
    dpObj <- ClusterParameterUpdate(dpObj)
    dpObj <- UpdateAlpha(dpObj)

    # Only update prior parameters for non-conjugate models when requested
    if (updatePrior && !inherits(dpObj$mixingDistribution, "conjugate")) {
      dpObj$mixingDistribution <- PriorParametersUpdate(dpObj$mixingDistribution,
                                                        dpObj$clusterParameters)
    }

    if (progressBar) {
      setTxtProgressBar(pb, i)
    }
  }

  dpObj$weights <- dpObj$pointsPerCluster / dpObj$n
  dpObj$alphaChain <- alphaChain
  dpObj$likelihoodChain <- likelihoodChain
  dpObj$weightsChain <- weightsChain
  dpObj$clusterParametersChain <- clusterParametersChain
  dpObj$priorParametersChain <- priorParametersChain
  dpObj$labelsChain <- labelsChain

  if (progressBar) {
    close(pb)
  }
  return(dpObj)
}

#' @export
Fit.conjugate <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive(), ...) {
  # Use C++ implementation if available and enabled
  if (using_cpp() && can_use_cpp(dpObj)) {
    return(Fit.dirichletprocess(dpObj, its, updatePrior, progressBar, ...))
  }

  # Otherwise use default R implementation
  return(Fit.default(dpObj, its, updatePrior, progressBar, ...))
}

#' @export
Fit.nonconjugate <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive(), ...) {
  # Use unified C++ implementation if available (including MVNormal2)
  if (using_cpp() && can_use_cpp(dpObj)) {
    return(Fit.dirichletprocess(dpObj, its, updatePrior, progressBar, ...))
  }

  # Otherwise use default R implementation
  return(Fit.default(dpObj, its, updatePrior, progressBar, ...))
}

#' @export
Fit.dirichletprocess <- function(dpObj, its, updatePrior = FALSE, progressBar = TRUE, ...) {
  # Validate inputs
  if (!inherits(dpObj, "dirichletprocess")) {
    stop("dpObj must be a dirichletprocess object")
  }
  if (its <= 0) {
    stop("Number of iterations must be positive")
  }

  # Extract additional parameters
  dots <- list(...)
  n_burn <- ifelse(is.null(dots$n_burn), 0, dots$n_burn)
  thin <- ifelse(is.null(dots$thin), 1, dots$thin)

  # Check if we should use C++ implementation
  use_cpp <- getOption("dirichletprocess.use_cpp", FALSE) && can_use_cpp(dpObj)

  if (use_cpp) {
    tryCatch({
      # Ensure dpObj has all required fields
      if (is.null(dpObj$data) || is.null(dpObj$alpha)) {
        stop("Invalid dirichletprocess object: missing data or alpha")
      }

      # Prepare parameters for C++
      mixing_params <- prepare_mixing_dist_params(dpObj)
      mcmc_params <- prepare_mcmc_params(dpObj, its, updatePrior, n_burn, thin)

      # Initialize cluster labels if not present
      if (is.null(dpObj$clusterLabels)) {
        dpObj$clusterLabels <- rep(1L, nrow(dpObj$data))
      }

      # Run C++ MCMC
      results <- run_mcmc_cpp(
        data = as.matrix(dpObj$data),
        mixing_dist_params = mixing_params,
        mcmc_params = mcmc_params
      )

      # Update dpObj with results
      if (!is.null(results$cluster_labels)) {
        # Get the final cluster labels and convert from 0-indexed to 1-indexed
        final_labels <- results$cluster_labels[[length(results$cluster_labels)]]
        dpObj$clusterLabels <- final_labels + 1
      }

      if (!is.null(results$alpha)) {
        # Get the final alpha value (handle both vector and list cases)
        alpha_chain <- results$alpha
        if (is.list(alpha_chain)) {
          dpObj$alpha <- as.numeric(alpha_chain[[length(alpha_chain)]])
        } else {
          dpObj$alpha <- as.numeric(tail(alpha_chain, 1))
        }
      }

      # Store chains (convert label chains from 0-indexed to 1-indexed)
      if (!is.null(results$labelsChain)) {
        dpObj$labelsChain <- lapply(results$labelsChain, function(labels) labels + 1)
      }
      dpObj$alphaChain <- results$alphaChain
      dpObj$likelihoodChain <- results$likelihoodChain

      # Extract final cluster parameters
      if (!is.null(results$theta_chain)) {
        final_params <- results$theta_chain[[length(results$theta_chain)]]
        
        # Convert parameter format for beta and beta2 distributions
        if (inherits(dpObj, "beta") || inherits(dpObj, "beta2")) {
          # C++ returns list(cluster1=c(mu1,nu1), cluster2=c(mu2,nu2), ...)
          # R expects list(mu=array(mu1,mu2,...), nu=array(nu1,nu2,...))
          n_clusters <- length(final_params)
          if (n_clusters > 0) {
            mu_vals <- sapply(final_params, function(x) x[1])
            nu_vals <- sapply(final_params, function(x) x[2])
            
            # Create arrays with proper dimensions for beta/beta2
            mu_array <- array(mu_vals, dim = c(1, 1, n_clusters))
            nu_array <- array(nu_vals, dim = c(1, 1, n_clusters))
            
            dpObj$clusterParameters <- list(mu = mu_array, nu = nu_array)
          }
        } else {
          dpObj$clusterParameters <- final_params
        }
      }

      # Update cluster counts
      unique_labels <- unique(dpObj$clusterLabels)
      dpObj$numberClusters <- length(unique_labels)
      dpObj$pointsPerCluster <- as.numeric(table(factor(dpObj$clusterLabels,
                                                        levels = seq_len(dpObj$numberClusters))))
      dpObj$weights <- dpObj$pointsPerCluster / dpObj$n

      # Store parameter chains
      if (inherits(dpObj, "beta") || inherits(dpObj, "beta2")) {
        # Convert parameter chain format for beta and beta2
        dpObj$clusterParametersChain <- lapply(results$theta_chain, function(iter_params) {
          n_clusters <- length(iter_params)
          if (n_clusters > 0) {
            mu_vals <- sapply(iter_params, function(x) x[1])
            nu_vals <- sapply(iter_params, function(x) x[2])
            
            # Create arrays with proper dimensions for beta/beta2
            mu_array <- array(mu_vals, dim = c(1, 1, n_clusters))
            nu_array <- array(nu_vals, dim = c(1, 1, n_clusters))
            
            list(mu = mu_array, nu = nu_array)
          } else {
            list(mu = array(dim = c(1, 1, 0)), nu = array(dim = c(1, 1, 0)))
          }
        })
      } else {
        dpObj$clusterParametersChain <- results$theta_chain
      }
      if (!is.null(results$cluster_labels)) {
        dpObj$weightsChain <- lapply(results$cluster_labels, function(labels) {
          # Convert 0-indexed to 1-indexed labels for weight calculation
          table(labels + 1) / length(labels)
        })
      }

      # Prior parameters chain if updated
      if (updatePrior && !is.null(results$prior_params_chain)) {
        dpObj$priorParametersChain <- results$prior_params_chain
        dpObj$mixingDistribution$priorParameters <-
          results$prior_params_chain[[length(results$prior_params_chain)]]
      }

      return(dpObj)

    }, error = function(e) {
      warning("C++ implementation failed: ", e$message,
              "\nFalling back to R implementation")
      return(Fit.default(dpObj, its, updatePrior, progressBar, ...))
    })
  }

  # Use R implementation
  return(Fit.default(dpObj, its, updatePrior, progressBar, ...))
}

#' @export
Fit.hierarchical <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive(), ...) {
  # Use C++ implementation if enabled and available
  if (using_cpp_hierarchical_samplers() && can_use_hierarchical_cpp(dpObj)) {
    return(Fit.hierarchical.cpp(dpObj, its, updatePrior, progressBar))
  }

  # Original R implementation
  if (progressBar) {
    pb <- txtProgressBar(min = 0, max = its, width = 50, char = "-", style = 3)
  }

  # Initialize storage arrays
  gammaValues <- numeric(its)
  gammaChain <- numeric(its)

  # Initialize alpha chains for each individual DP
  for (j in seq_along(dpObj$indDP)) {
    dpObj$indDP[[j]]$alphaChain <- numeric(its)
    dpObj$indDP[[j]]$likelihoodChain <- numeric(its)
    dpObj$indDP[[j]]$weightsChain <- vector("list", length = its)
    dpObj$indDP[[j]]$clusterParametersChain <- vector("list", length = its)
    dpObj$indDP[[j]]$labelsChain <- vector("list", length = its)
  }

  # Initialize global parameter storage
  globalParametersChain <- vector("list", length = its)
  globalStickChain <- vector("list", length = its)

  for (i in seq_len(its)) {
    # Update cluster components for each individual DP
    dpObj <- ClusterComponentUpdate(dpObj)

    # Update cluster parameters for each individual DP
    dpObj <- ClusterParameterUpdate(dpObj)

    # Update alpha for each individual DP
    dpObj <- UpdateAlpha(dpObj)

    # Update global parameters using all data
    dpObj <- GlobalParameterUpdate(dpObj)

    # Update G0 (the base distribution)
    dpObj <- UpdateG0(dpObj)

    # Update gamma (concentration parameter for G0)
    dpObj <- UpdateGamma(dpObj)

    # Store values for this iteration
    gammaValues[i] <- dpObj$gamma
    gammaChain[i] <- dpObj$gamma
    globalParametersChain[[i]] <- dpObj$globalParameters
    globalStickChain[[i]] <- dpObj$globalStick

    # Store individual DP values
    for (j in seq_along(dpObj$indDP)) {
      # Store alpha
      dpObj$indDP[[j]]$alphaChain[i] <- dpObj$indDP[[j]]$alpha

      # Calculate and store likelihood
      if (!is.null(dpObj$indDP[[j]]$data) && !is.null(dpObj$indDP[[j]]$clusterLabels)) {
        dpObj$indDP[[j]]$likelihoodChain[i] <- sum(log(LikelihoodDP(dpObj$indDP[[j]])))
      }

      # Store weights
      dpObj$indDP[[j]]$weightsChain[[i]] <- dpObj$indDP[[j]]$pointsPerCluster / dpObj$indDP[[j]]$n

      # Store cluster parameters
      dpObj$indDP[[j]]$clusterParametersChain[[i]] <- dpObj$indDP[[j]]$clusterParameters

      # Store labels
      dpObj$indDP[[j]]$labelsChain[[i]] <- dpObj$indDP[[j]]$clusterLabels

      # Update weights
      dpObj$indDP[[j]]$weights <- dpObj$indDP[[j]]$pointsPerCluster / dpObj$indDP[[j]]$n
    }

    # Update prior parameters if requested
    if (updatePrior) {
      # Get unique cluster parameters across all DPs
      allClusterParams <- list()
      for (j in seq_along(dpObj$indDP)) {
        if (!is.null(dpObj$indDP[[j]]$clusterParameters)) {
          allClusterParams <- c(allClusterParams,
                                list(dpObj$indDP[[j]]$clusterParameters))
        }
      }

      # Find unique parameters
      if (length(allClusterParams) > 0) {
        uniqueParams <- unique(unlist(lapply(allClusterParams, function(x) {
          if (is.list(x)) x[[1]] else x
        }), recursive = FALSE))

        clustParamLen <- length(uniqueParams)

        if (clustParamLen > 0) {
          # Extract global parameters up to the number of unique clusters
          clustParam <- lapply(dpObj$globalParameters, function(x) {
            if (is.array(x) && length(dim(x)) >= 3) {
              x[, , 1:min(clustParamLen, dim(x)[3]), drop = FALSE]
            } else {
              x
            }
          })

          # Update prior parameters using the first DP's mixing distribution
          tempMD <- PriorParametersUpdate(dpObj$indDP[[1]]$mixingDistribution, clustParam)

          # Apply updated prior parameters to all individual DPs
          for (j in seq_along(dpObj$indDP)) {
            dpObj$indDP[[j]]$mixingDistribution$priorParameters <- tempMD$priorParameters
          }
        }
      }
    }

    if (progressBar) {
      setTxtProgressBar(pb, i)
    }
  }

  # Store all chains in the dpObj
  dpObj$gammaValues <- gammaValues
  dpObj$gammaChain <- gammaChain
  dpObj$globalParametersChain <- globalParametersChain
  dpObj$globalStickChain <- globalStickChain

  # Ensure each individual DP has the correct numberClusters as a scalar
  for (j in seq_along(dpObj$indDP)) {
    if (!is.null(dpObj$indDP[[j]]$clusterLabels)) {
      dpObj$indDP[[j]]$numberClusters <- length(unique(dpObj$indDP[[j]]$clusterLabels))
    }
  }

  if (progressBar) {
    close(pb)
  }

  return(dpObj)
}

#' @export
Fit.hierarchical.cpp <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive(), ...) {
  if (!can_use_hierarchical_cpp(dpObj)) {
    stop("C++ implementation not available for this hierarchical DP type")
  }

  # Use the appropriate C++ implementation based on distribution type
  if (all(sapply(dpObj$indDP, function(x) inherits(x, "beta")))) {
    # Use hierarchical Beta C++ implementation
    result <- run_hierarchical_mcmc_cpp(
      dpObj,
      n_iter = its,
      n_burn = 0,  # No burn-in for regular Fit
      thin = 1,
      update_prior = updatePrior,
      progress_bar = progressBar
    )
  } else if (all(sapply(dpObj$indDP, function(x) inherits(x, "mvnormal2")))) {
    # Use hierarchical MVNormal2 C++ implementation (disable progress bar to avoid R implementation)
    result <- hierarchical_mvnormal2_fit_cpp(
      dpObj,
      iterations = its,
      updatePrior = updatePrior,
      progressBar = FALSE  # Disable to ensure C++ is used
    )
  } else {
    stop("Mixed distribution types not supported in hierarchical C++ implementation")
  }

  # The result from run_hierarchical_mcmc_cpp should already have the updated dpObj
  # Ensure all fields are properly set

  # Make sure numberClusters is scalar for each individual DP
  for (j in seq_along(result$indDP)) {
    if (!is.null(result$indDP[[j]]$clusterLabels)) {
      result$indDP[[j]]$numberClusters <- as.integer(length(unique(result$indDP[[j]]$clusterLabels)))
    }

    # Ensure weights are calculated
    if (!is.null(result$indDP[[j]]$pointsPerCluster) && !is.null(result$indDP[[j]]$n)) {
      result$indDP[[j]]$weights <- result$indDP[[j]]$pointsPerCluster / result$indDP[[j]]$n
    }
  }

  # Ensure gamma is set to the last value if we have samples
  if (!is.null(result$gammaValues) && length(result$gammaValues) > 0) {
    result$gamma <- result$gammaValues[length(result$gammaValues)]
  }

  # Set gammaChain as alias for gammaValues for compatibility
  if (!is.null(result$gammaValues)) {
    result$gammaChain <- result$gammaValues
  }

  return(result)
}

#' @export
Fit.markov <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive(), ...) {
  # Similar pattern - check for C++ then fall back to R
  if (using_cpp() && exists("_dirichletprocess_markov_dp_fit_cpp")) {
    return(Fit.markov.cpp(dpObj, its, updatePrior, progressBar))
  }

  # R implementation would go here
  return(dpObj)
}
