#' C++ Implementation Wrappers for Markov DP (HMM)
#'
#' These functions provide access to the C++ implementations of the
#' Markov Dirichlet process (Hidden Markov Model) algorithms.
#'
#' @name cpp_markov_wrappers
#' @keywords internal
NULL

#' @rdname cpp_markov_wrappers
#' @param dp Markov Dirichlet process object
#' @param its Number of iterations
#' @param progressBar Whether to show progress bar
#' @export
Fit.markov.cpp <- function(dp, its, updatePrior = FALSE, progressBar = TRUE, ...) {
  if (!inherits(dp, "markov")) {
    stop("This C++ implementation is only for Markov Dirichlet processes")
  }

  # Convert 1-indexed R states to 0-indexed C++ states
  dp$states <- dp$states - 1

  # Call C++ implementation
  result <- markov_dp_fit_cpp(dp, its, updatePrior, progressBar)

  # Convert back to 1-indexed
  result$states <- result$states + 1
  if (!is.null(result$statesChain)) {
    result$statesChain <- lapply(result$statesChain, function(x) x + 1)
  }

  return(result)
}

#' @rdname cpp_markov_wrappers
#' @export
UpdateStates.cpp <- function(dp) {
  if (!inherits(dp, "markov")) {
    stop("This C++ implementation is only for Markov Dirichlet processes")
  }

  # Convert states
  dp$states <- dp$states - 1

  # Call C++ implementation
  result <- markov_dp_update_states_cpp(dp)

  # Convert back
  result$states <- result$states + 1

  # Update the params list based on new states
  new_states <- list(result$states, result$params)

  return(new_states)
}

#' @rdname cpp_markov_wrappers
#' @export
UpdateAlphaBeta.cpp <- function(dp) {
  if (!inherits(dp, "markov")) {
    stop("This C++ implementation is only for Markov Dirichlet processes")
  }

  # Convert states
  dp$states <- dp$states - 1

  # Call C++ implementation
  result <- markov_dp_update_alpha_beta_cpp(dp)

  # Convert back
  result$states <- result$states + 1

  # Extract alpha and beta
  dp$alpha <- result$alpha
  dp$beta <- result$beta

  return(dp)
}

#' @rdname cpp_markov_wrappers
#' @export
param_update.cpp <- function(dp) {
  if (!inherits(dp, "markov")) {
    stop("This C++ implementation is only for Markov Dirichlet processes")
  }

  # Convert states
  dp$states <- dp$states - 1

  # Call C++ implementation
  result <- markov_dp_param_update_cpp(dp)

  # Convert back
  result$states <- result$states + 1

  # Update dp with new parameters
  dp$uniqueParams <- result$uniqueParams
  dp$params <- result$params

  return(dp)
}

#' Enable C++ implementations for Markov DP samplers
#'
#' This function enables the use of C++ implementations for the Markov
#' DP (HMM) sampling algorithms when available.
#'
#' @param use_cpp Logical indicating whether to use C++ implementations
#' @export
enable_cpp_markov_samplers <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp_markov = use_cpp)

  if (use_cpp) {
    message("C++ samplers enabled for Markov Dirichlet processes (HMM)")
  } else {
    message("Using R implementations for Markov samplers")
  }

  invisible(use_cpp)
}

#' Check if C++ Markov samplers are enabled
#'
#' @return Logical indicating if C++ Markov samplers are enabled
#' @export
using_cpp_markov_samplers <- function() {
  getOption("dirichletprocess.use_cpp_markov", FALSE)
}
