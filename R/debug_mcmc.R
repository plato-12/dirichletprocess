# R/debug_mcmc.R
#' Debug MCMC C++ implementation
#' @export
debug_mcmc_cpp <- function(data, n_iter = 10, verbose = TRUE) {
  if (!exists("_dirichletprocess_run_mcmc_cpp")) {
    stop("C++ implementation not compiled")
  }

  data_matrix <- matrix(data, ncol = 1)

  mcmc_params <- list(
    n_iter = n_iter,
    n_burn = 0,
    thin = 1,
    update_concentration = TRUE,
    alpha = 1.0
  )

  dist_params <- list(
    type = "gaussian",
    mu0 = mean(data),
    kappa0 = 0.01,
    alpha0 = 2.0,
    beta0 = var(data)
  )

  result <- run_mcmc_cpp(data_matrix, dist_params, mcmc_params)

  if (verbose) {
    cat("Debug MCMC Results:\n")
    cat("Final clusters:", result$final_n_clusters, "\n")
    cat("Cluster evolution:", result$n_clusters, "\n")
    cat("Alpha evolution:", round(result$alpha, 3), "\n")
  }

  return(result)
}
