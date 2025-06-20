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
      mcmc_params <- prepare_mcmc_params(dpObj, its, updatePrior, n_burn, thin)

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
        # Get the final cluster labels
        dpObj$clusterLabels <- results$cluster_labels[[length(results$cluster_labels)]]
      }

      if (!is.null(results$alpha)) {
        # Get the final alpha value
        dpObj$alpha <- tail(results$alpha, 1)[[1]][1]
      }

      # Store chains
      dpObj$labelsChain <- results$labelsChain
      dpObj$alphaChain <- results$alphaChain
      dpObj$likelihoodChain <- results$likelihoodChain

      # Update cluster parameters
      if (!is.null(results$theta) && length(results$theta) > 0) {
        final_theta <- results$theta[[length(results$theta)]]
        dpObj$clusterParameters <- final_theta
      }

      # Update number of clusters
      dpObj$numberClusters <- length(unique(dpObj$clusterLabels))

      # Update points per cluster
      dpObj$pointsPerCluster <- as.numeric(table(dpObj$clusterLabels))

      return(dpObj)

    }, error = function(e) {
      if (progressBar) {
        message("C++ implementation failed, falling back to R: ", e$message)
      }
      use_cpp <- FALSE
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
  # Check if C++ implementation exists
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    return(FALSE)
  }

  # Add mvnormal to supported types
  supported_types <- c("normal_inverse_gamma", "normal", "gaussian", "beta", "mvnormal")

  # Check if the mixing distribution inherits from any supported type
  return(inherits(dp_obj$mixingDistribution, supported_types))
}

#' @export
Fit.conjugate <- function(dpObj, its = 1000, updatePrior = TRUE, progressBar = TRUE, ...) {

  if (its <= 0) {
    stop("Number of iterations must be positive")
  }

  # Check if we should use C++ implementation
  use_cpp <- getOption("dirichletprocess.use_cpp", FALSE) && can_use_cpp(dpObj)

  if (use_cpp && inherits(dpObj$mixingDistribution, "mvnormal")) {
    # MVNormal has specialized C++ functions
    return(Fit_mvnormal_cpp(dpObj, its, updatePrior, progressBar))
  }

  # Otherwise use default implementation
  return(Fit.default(dpObj, its, updatePrior, progressBar, ...))
}

#' Specialized Fit function for MVNormal with C++
#' @keywords internal
Fit_mvnormal_cpp <- function(dpObj, its, updatePrior, progressBar) {

  # Check if MVNormal C++ functions are available
  if (!exists("conjugate_mvnormal_cluster_component_update_cpp") ||
      !exists("conjugate_mvnormal_cluster_parameter_update_cpp")) {
    message("MVNormal C++ functions not available, falling back to R")
    return(Fit.default(dpObj, its, updatePrior, progressBar))
  }

  # Initialize chains
  dpObj$alphaChain <- numeric(its)
  dpObj$clusterLabelChain <- vector("list", its)
  dpObj$likelihoodChain <- numeric(its)

  if (progressBar) {
    pb <- txtProgressBar(min = 0, max = its, style = 3)
  }

  for (i in seq_len(its)) {
    # Update cluster assignments using C++
    dpObj <- ClusterComponentUpdate.mvnormal.cpp(dpObj)

    # Update cluster parameters using C++
    dpObj <- ClusterParameterUpdate.mvnormal.cpp(dpObj)

    # Update alpha (using R for now)
    if (updatePrior) {
      dpObj <- UpdateAlpha(dpObj)
    }

    # Store iteration results
    dpObj$alphaChain[i] <- dpObj$alpha
    dpObj$clusterLabelChain[[i]] <- dpObj$clusterLabels
    dpObj$likelihoodChain[i] <- sum(log(dpObj$pointsPerCluster))

    if (progressBar) {
      setTxtProgressBar(pb, i)
    }
  }

  if (progressBar) {
    close(pb)
  }

  dpObj$iterations <- its
  return(dpObj)
}
