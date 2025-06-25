#' Run Hierarchical Beta MCMC using C++ implementation
#'
#' @param dp_list Hierarchical DP object (not a list of DirichletProcessBeta objects)
#' @param n_iter Number of MCMC iterations
#' @param n_burn Number of burn-in iterations
#' @param thin Thinning parameter
#' @param update_prior Whether to update prior parameters
#' @param progress_bar Show progress bar
#'
#' @return Updated hierarchical DP object
#' @export
run_hierarchical_mcmc_cpp <- function(dp_list, n_iter = 1000, n_burn = 100,
                                      thin = 1, update_prior = FALSE,
                                      progress_bar = TRUE) {

  # Validate inputs - dp_list should be a hierarchical DP object
  if (!inherits(dp_list, "hierarchical")) {
    stop("dp_list must be a hierarchical Dirichlet process object")
  }

  if (!all(sapply(dp_list$indDP, function(x) inherits(x, "beta")))) {
    stop("All individual DPs must be Beta type")
  }

  # Extract datasets
  datasets <- lapply(dp_list$indDP, function(dp) as.matrix(dp$data))

  # Prepare mixing distribution parameters
  first_dp <- dp_list$indDP[[1]]
  mixing_params <- list(
    type = "hierarchical_beta",
    alpha0 = first_dp$mixingDistribution$priorParameters[1],
    beta0 = first_dp$mixingDistribution$priorParameters[2],
    maxT = first_dp$mixingDistribution$maxT,
    gamma_prior_shape = dp_list$gammaPriors[1],
    gamma_prior_rate = dp_list$gammaPriors[2]
  )

  # MCMC parameters - use the first alpha value as default
  # The C++ code expects a single 'alpha' parameter
  mcmc_params <- list(
    n_iter = as.integer(n_iter),
    n_burn = as.integer(n_burn),
    thin = as.integer(thin),
    update_prior = update_prior,
    update_concentration = TRUE,
    m_auxiliary = 3L,  # For Algorithm 8
    alpha = as.numeric(dp_list$indDP[[1]]$alpha)  # Use first DP's alpha as default
  )

  # Call C++ implementation
  result <- .Call("_dirichletprocess_run_hierarchical_mcmc_cpp",
                  datasets, mixing_params, mcmc_params,
                  PACKAGE = "dirichletprocess")

  # Update dp_list with results
  for (i in seq_along(result$individual_dps)) {
    dp_list$indDP[[i]] <- update_dp_from_mcmc(
      dp_list$indDP[[i]],
      result$individual_dps[[i]]
    )
  }

  dp_list$globalParameters <- result$global_parameters
  dp_list$globalStick <- result$global_weights
  dp_list$gammaValues <- result$gamma_samples

  # Update gamma if it was sampled
  if (length(result$gamma_samples) > 0) {
    dp_list$gamma <- result$gamma_samples[length(result$gamma_samples)]
  }

  return(dp_list)
}

#' Check if hierarchical C++ implementation is available
#'
#' @param dp_list Hierarchical DP object
#' @return Logical indicating availability
#' @export
can_use_hierarchical_cpp <- function(dp_list) {
  if (!exists("_dirichletprocess_run_hierarchical_mcmc_cpp")) {
    return(FALSE)
  }

  # Check if it's a hierarchical object
  if (!inherits(dp_list, "hierarchical")) {
    return(FALSE)
  }

  # Check if all individual DPs are supported
  all_beta <- all(sapply(dp_list$indDP, function(x) inherits(x, "beta")))

  return(all_beta)
}

#' Update DP object from MCMC results
#' @param dp Original DP object
#' @param mcmc_result MCMC results from C++
#' @return Updated DP object
#' @keywords internal
update_dp_from_mcmc <- function(dp, mcmc_result) {
  if (!is.null(mcmc_result$cluster_labels)) {
    dp$clusterLabels <- mcmc_result$cluster_labels
  }

  if (!is.null(mcmc_result$cluster_params)) {
    dp$clusterParameters <- mcmc_result$cluster_params
  }

  if (!is.null(mcmc_result$n_clusters)) {
    dp$numberClusters <- mcmc_result$n_clusters
  }

  if (!is.null(mcmc_result$alpha)) {
    dp$alpha <- mcmc_result$alpha
  }

  if (!is.null(mcmc_result$weights)) {
    dp$weights <- mcmc_result$weights
  }

  return(dp)
}
