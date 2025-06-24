#' Run Hierarchical Beta MCMC using C++ implementation
#'
#' @param dp_list List of DirichletProcessBeta objects
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

  # Validate inputs
  if (!all(sapply(dp_list$indDP, function(x) inherits(x, "beta")))) {
    stop("All individual DPs must be Beta type")
  }

  # Extract datasets
  datasets <- lapply(dp_list$indDP, function(dp) dp$data)

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

  # MCMC parameters
  mcmc_params <- list(
    n_iter = n_iter,
    n_burn = n_burn,
    thin = thin,
    update_prior = update_prior,
    update_concentration = TRUE,
    m_auxiliary = 3  # For Algorithm 8
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

  # Check if all individual DPs are supported
  all_beta <- all(sapply(dp_list$indDP, function(x) inherits(x, "beta")))

  return(all_beta)
}
