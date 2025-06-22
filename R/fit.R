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
  # For nonconjugate, check if C++ implementation is available
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
        # Get the final cluster labels
        dpObj$clusterLabels <- results$cluster_labels[[length(results$cluster_labels)]]
      }

      if (!is.null(results$alpha)) {
        # Get the final alpha value
        dpObj$alpha <- tail(results$alpha, 1)
      }

      # Store chains
      dpObj$labelsChain <- results$labelsChain
      dpObj$alphaChain <- results$alphaChain
      dpObj$likelihoodChain <- results$likelihoodChain

      # Extract final cluster parameters
      if (!is.null(results$cluster_params)) {
        final_params <- results$cluster_params[[length(results$cluster_params)]]
        dpObj$clusterParameters <- final_params
      }

      # Update cluster counts
      unique_labels <- unique(dpObj$clusterLabels)
      dpObj$numberClusters <- length(unique_labels)
      dpObj$pointsPerCluster <- as.numeric(table(factor(dpObj$clusterLabels,
                                                        levels = seq_len(dpObj$numberClusters))))
      dpObj$weights <- dpObj$pointsPerCluster / dpObj$n

      # Store parameter chains
      dpObj$clusterParametersChain <- results$cluster_params
      dpObj$weightsChain <- lapply(results$cluster_labels, function(labels) {
        table(labels) / length(labels)
      })

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
Fit.markov <- function(dpObj, its = 1000, updatePrior = FALSE, progressBar = interactive(), ...) {
  # Similar pattern - check for C++ then fall back to R
  if (using_cpp() && exists("_dirichletprocess_markov_dp_fit_cpp")) {
    return(Fit.markov.cpp(dpObj, its, updatePrior, progressBar))
  }

  # R implementation would go here
  return(dpObj)
}
