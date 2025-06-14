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
Fit.dirichletprocess <- function(dp_obj, n_iter, n_burn = 0, thin = 1,
                                 update_concentration = TRUE, ...) {

  if (getOption("dirichletprocess.use_cpp", FALSE) &&
      can_use_cpp(dp_obj)) {

    # Prepare parameters for C++
    mixing_params <- prepare_mixing_dist_params(dp_obj)
    mcmc_params <- list(
      n_iter = n_iter,
      n_burn = n_burn,
      thin = thin,
      update_concentration = update_concentration,
      alpha = dp_obj$alpha
    )

    # Run C++ MCMC
    results <- run_mcmc_cpp(dp_obj$data, mixing_params, mcmc_params)

    # Update dp_obj with results
    dp_obj$cluster_labels <- results$cluster_labels
    dp_obj$alpha <- results$alpha
    dp_obj$theta <- results$theta
    dp_obj$n_clusters <- results$n_clusters

  } else {
    # FIXED: Use correct R implementation that doesn't confuse update_concentration with updatePrior
    dp_obj <- fit_r_implementation(dp_obj, n_iter, n_burn, thin,
                                   update_concentration, ...)
  }

  return(dp_obj)
}

#' R implementation function for backward compatibility and testing
#' @keywords internal
fit_r_implementation <- function(dp_obj, n_iter, n_burn = 0, thin = 1,
                                 update_concentration = TRUE, ...) {

  # FIXED: Extract and handle parameters properly to avoid conflicts
  dots <- list(...)

  # Remove updatePrior from dots if it exists to avoid conflict
  if ("updatePrior" %in% names(dots)) {
    updatePrior <- dots$updatePrior
    dots$updatePrior <- NULL
  } else {
    # For conjugate models, we typically don't update prior parameters
    updatePrior <- FALSE
  }

  # If this is a non-conjugate model, we might want to update prior parameters
  if (!inherits(dp_obj$mixingDistribution, "conjugate")) {
    # For now, keeping it FALSE unless explicitly requested
    if (!"updatePrior" %in% names(list(...))) {
      updatePrior <- FALSE
    }
  }

  # Call the standard R implementation with correct parameters
  # Pass remaining arguments from dots
  do.call(Fit.default, c(list(dpObj = dp_obj, its = n_iter, updatePrior = updatePrior), dots))
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
