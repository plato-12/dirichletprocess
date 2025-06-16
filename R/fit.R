#' Fit the Dirichlet process object
#'
#' Using Neal's algorithm 4 or 8 depending on conjugacy the sampling procedure for a Dirichlet process is carried out.
#' Lists of both cluster parameters, weights and the sampled concentration values are included in the fitted \code{dpObj}.
#' When \code{update_prior} is set to \code{TRUE} the parameters of the base measure are also updated.
#'
#' @param dpObj Initialised Dirichlet Process object
#' @param its Number of iterations to use
#' @param updatePrior Logical flag, defaults to \code{FALSE}. Set whether the parameters of the base measure are updated.
#' @param progressBar Logical flag indicating whether to display a progress bar.
#' @return A Dirichlet Process object with the fitted cluster parameters and labels.
#'
#' @references Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture models. Journal of computational and graphical statistics, 9(2), 249-265.
#'
#' @export
Fit <- function(dpObj, its, updatePrior = FALSE, progressBar = TRUE) UseMethod("Fit", dpObj)

#' @export
Fit.default <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive()) {

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

    # FIXED: Only update prior parameters for non-conjugate models and when explicitly requested
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
Fit.hierarchical <- function(dpObj, its, updatePrior = FALSE, progressBar = interactive()) {
  # Use C++ implementation if enabled and available
  if (using_cpp_hierarchical_samplers() && all(sapply(dpObj$indDP, function(x) inherits(x, "beta")))) {
    return(Fit.hierarchical.cpp(dpObj, its, updatePrior, progressBar))
  }

  # Original R implementation
  if (progressBar) {
    pb <- txtProgressBar(min = 0, max = its, width = 50, char = "-", style = 3)
  }

  gammaValues <- numeric(its)

  for (i in seq_len(its)) {

    dpObj <- ClusterComponentUpdate(dpObj)
    dpObj <- UpdateAlpha(dpObj)
    dpObj <- GlobalParameterUpdate(dpObj)
    dpObj <- UpdateG0(dpObj)
    dpObj <- UpdateGamma(dpObj)

    if (updatePrior) {

      clustParamLen <- length(unique(lapply(dpObj$indDP, function(x) x$clusterParameters[[1]])))

      clustParam <- lapply(dpObj$globalParameters, function(x) x[, , 1:clustParamLen, drop = FALSE])

      tempMD <- PriorParametersUpdate(dpObj$indDP[[1]]$mixingDistribution, clustParam)

      for (j in seq_along(dpObj$indDP)) {
        dpObj$indDP[[j]]$mixingDistribution$priorParameters <- tempMD$priorParameters
      }
    }

    if (progressBar) {
      setTxtProgressBar(pb, i)
    }

    gammaValues[i] <- dpObj$gamma
  }
  dpObj$gammaValues <- gammaValues
  if (progressBar) {
    close(pb)
  }
  return(dpObj)
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

      mcmc_params <- list(
        n_iter = as.integer(its),
        n_burn = as.integer(n_burn),
        thin = as.integer(thin),
        update_concentration = as.logical(updatePrior),
        alpha = as.numeric(dpObj$alpha),
        progressBar = as.logical(progressBar)
      )

      # Initialize cluster labels if not present
      if (is.null(dpObj$clusterLabels)) {
        dpObj$clusterLabels <- rep(1L, length(dpObj$data))
      }

      # Run C++ MCMC
      results <- run_mcmc_cpp(
        data = as.matrix(dpObj$data),
        mixing_dist_params = mixing_params,
        mcmc_params = mcmc_params
      )

      # Update dpObj with results
      if (!is.null(results$cluster_labels)) {
        # Get the final cluster labels (last iteration)
        final_labels <- if (is.matrix(results$cluster_labels)) {
          results$cluster_labels[nrow(results$cluster_labels), ]
        } else {
          results$cluster_labels
        }
        dpObj$clusterLabels <- as.integer(final_labels)
      }

      if (!is.null(results$alpha)) {
        # Get the final alpha value
        dpObj$alpha <- if (length(results$alpha) > 1) {
          tail(results$alpha, 1)
        } else {
          results$alpha
        }
      }

      # Update cluster parameters and related fields
      if (!is.null(results$theta)) {
        dpObj$clusterParameters <- results$theta
      }

      # Calculate cluster statistics
      if (!is.null(dpObj$clusterLabels)) {
        dpObj$pointsPerCluster <- as.integer(table(dpObj$clusterLabels))
        dpObj$numberClusters <- length(unique(dpObj$clusterLabels))
        dpObj$weights <- dpObj$pointsPerCluster / length(dpObj$data)
      }

      # Store chains if available
      if (!is.null(results$alpha_chain)) {
        dpObj$alphaChain <- results$alpha_chain
      } else if (!is.null(results$alpha) && length(results$alpha) > 1) {
        dpObj$alphaChain <- results$alpha
      }

      if (!is.null(results$likelihood_chain)) {
        dpObj$likelihoodChain <- results$likelihood_chain
      }

      if (!is.null(results$labels_chain)) {
        dpObj$labelsChain <- results$labels_chain
      }

      if (!is.null(results$theta_chain)) {
        dpObj$clusterParametersChain <- results$theta_chain
      }

    }, error = function(e) {
      warning("C++ implementation failed: ", e$message, ". Falling back to R implementation.")
      dpObj <- Fit.default(dpObj = dpObj, its = its,
                           updatePrior = updatePrior,
                           progressBar = progressBar)
    })
  } else {
    # Use R implementation
    dpObj <- Fit.default(dpObj = dpObj, its = its,
                         updatePrior = updatePrior,
                         progressBar = progressBar)
  }

  # Set the iterations count for the print method
  dpObj$iterations <- its

  return(dpObj)
}

#' R implementation function for backward compatibility and testing
#' @keywords internal
fit_r_implementation <- function(dp_obj, n_iter, n_burn = 0, thin = 1,
                                 update_concentration = TRUE, progressBar = TRUE, ...) {
  # FIXED: Extract and handle parameters properly to avoid conflicts
  dots <- list(...)
  # Remove updatePrior from dots if it exists to avoid conflict
  if ("updatePrior" %in% names(dots)) {
    updatePrior <- dots$updatePrior
    dots$updatePrior <- NULL
  } else {
    # For conjugate models, we typically don't update prior parameters
    updatePrior <- update_concentration
  }

  # If this is a non-conjugate model, we might want to update prior parameters
  if (!inherits(dp_obj$mixingDistribution, "conjugate")) {
    # For now, keeping it FALSE unless explicitly requested
    if (!"updatePrior" %in% names(list(...))) {
      updatePrior <- FALSE
    }
  }

  # Call the standard R implementation with correct parameters
  # Pass progressBar and remaining arguments from dots
  do.call(Fit.default, c(list(dpObj = dp_obj, its = n_iter,
                              updatePrior = updatePrior,
                              progressBar = progressBar), dots))
}

#' Check if C++ implementation is available for this model
#' @keywords internal
can_use_cpp <- function(dp_obj) {
  # For now, only return TRUE if we actually have the C++ implementation
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    return(FALSE)
  }

  supported_types <- c("normal_inverse_gamma", "normal")
  inherits(dp_obj$mixingDistribution, supported_types)
}
