# Utility functions specific to Beta distribution benchmarking

#' Run Beta scaling analysis
#'
#' @param max_n Maximum number of observations
#' @param n_points Number of data points to test
#' @return Data frame with scaling results
analyze_beta_scaling <- function(max_n = 5000, n_points = 10) {

  ns <- round(exp(seq(log(50), log(max_n), length.out = n_points)))

  results <- data.frame(
    n = integer(),
    implementation = character(),
    time = numeric(),
    time_per_obs = numeric(),
    time_per_iter = numeric()
  )

  n_iter <- 50  # Fixed iterations for scaling test

  for (n in ns) {
    cat(sprintf("Testing n = %d\n", n))

    # Generate data
    data <- generate_beta_mixture(n, list(c(2, 8), c(8, 2)))

    # Test R implementation
    set_use_cpp(FALSE)
    time_r <- system.time({
      dp <- DirichletProcessBeta(data, verbose = FALSE)
      dp <- Fit(dp, n_iter, progressBar = FALSE)
    })["elapsed"]

    results <- rbind(results, data.frame(
      n = n,
      implementation = "R",
      time = time_r,
      time_per_obs = time_r / n,
      time_per_iter = time_r / n_iter
    ))

    # Test C++ if available
    if (can_use_cpp(DirichletProcessBeta(data[1:10]))) {
      set_use_cpp(TRUE)
      time_cpp <- system.time({
        dp <- DirichletProcessBeta(data, verbose = FALSE)
        dp <- Fit(dp, n_iter, progressBar = FALSE)
      })["elapsed"]

      results <- rbind(results, data.frame(
        n = n,
        implementation = "C++",
        time = time_cpp,
        time_per_obs = time_cpp / n,
        time_per_iter = time_cpp / n_iter
      ))
    }
  }

  return(results)
}

#' Analyze Beta MH acceptance rates
#'
#' @param n_obs Number of observations
#' @param n_iter Number of iterations
#' @return List with acceptance rate analysis
analyze_beta_mh_acceptance <- function(n_obs = 100, n_iter = 1000) {

  # Generate challenging data (near boundaries)
  data <- c(rbeta(n_obs/2, 0.5, 10), rbeta(n_obs/2, 10, 0.5))

  # Create custom DP object to track acceptances
  dp <- DirichletProcessBeta(data, verbose = FALSE)

  # Storage for tracking
  acceptances <- numeric(n_iter)
  param_trace <- matrix(NA, nrow = n_iter, ncol = 2 * dp$numberClusters)

  # Run MCMC with tracking
  for (i in 1:n_iter) {
    old_params <- dp$clusterParameters
    dp <- ClusterParameterUpdate(dp)
    new_params <- dp$clusterParameters

    # Check if parameters changed (indicates acceptance)
    accepted <- !identical(old_params, new_params)
    acceptances[i] <- as.numeric(accepted)

    # Store parameters
    for (k in 1:dp$numberClusters) {
      param_trace[i, (2*k-1):(2*k)] <- c(
        dp$clusterParameters[[1]][k],
        dp$clusterParameters[[2]][k]
      )
    }
  }

  # Calculate acceptance rates
  overall_rate <- mean(acceptances)
  window_rates <- sapply(seq(1, n_iter-99, by = 10), function(i) {
    mean(acceptances[i:(i+99)])
  })

  return(list(
    overall_rate = overall_rate,
    window_rates = window_rates,
    param_trace = param_trace,
    acceptances = acceptances
  ))
}

#' Compare Beta implementations on edge cases
#'
#' @return Data frame with edge case results
benchmark_beta_edge_cases <- function() {

  edge_cases <- list(
    "Near zero" = rep(0.001, 50),
    "Near one" = rep(0.999, 50),
    "Bimodal extreme" = c(rep(0.001, 25), rep(0.999, 25)),
    "Uniform" = runif(50),
    "High precision" = rbeta(50, 100, 100),
    "Low precision" = rbeta(50, 0.1, 0.1)
  )

  results <- list()

  for (case_name in names(edge_cases)) {
    cat(sprintf("Testing edge case: %s\n", case_name))

    data <- edge_cases[[case_name]]

    # R implementation
    set_use_cpp(FALSE)
    time_r <- system.time({
      dp_r <- DirichletProcessBeta(data, verbose = FALSE)
      dp_r <- Fit(dp_r, 50, progressBar = FALSE)
    })["elapsed"]

    success_r <- !is.null(dp_r) && all(is.finite(dp_r$clusterParameters[[1]]))

    # C++ implementation
    if (can_use_cpp(DirichletProcessBeta(data[1:min(10, length(data))]))) {
      set_use_cpp(TRUE)
      time_cpp <- system.time({
        dp_cpp <- DirichletProcessBeta(data, verbose = FALSE)
        dp_cpp <- Fit(dp_cpp, 50, progressBar = FALSE)
      })["elapsed"]

      success_cpp <- !is.null(dp_cpp) && all(is.finite(dp_cpp$clusterParameters[[1]]))
    } else {
      time_cpp <- NA
      success_cpp <- NA
    }

    results[[case_name]] <- data.frame(
      case = case_name,
      implementation = c("R", "C++"),
      time = c(time_r, time_cpp),
      success = c(success_r, success_cpp)
    )
  }

  return(bind_rows(results))
}
