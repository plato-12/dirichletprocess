#' Create a generic Dirichlet process hidden Markov Model
#'
#' Create a hidden Markov model where the data is believed to be generated from the mixing object distribution.
#'
#' @param x Data to be modelled
#' @param mdobj Mixing distribution object
#' @param alpha Alpha parameter
#' @param beta Beta parameter
#' @param cpp Logical. Use C++ implementation if TRUE, R implementation if FALSE. Default is FALSE.
#' @export
DirichletHMMCreate <- function(x, mdobj, alpha, beta, cpp = FALSE){

  if(is.vector(x)){
    x <- matrix(x, ncol=1)
  }

  states <- seq_len(nrow(x))

  # Draw initial parameters for each state
  uniqueParams <- PriorDraw(mdobj, length(unique(states)))

  # Create params list with proper structure
  params <- vector("list", length(states))

  for (i in seq_along(states)) {
    state_idx <- states[i]
    params[[i]] <- list()

    # Extract parameters for this state with proper structure
    for (j in seq_along(uniqueParams)) {
      param_array <- uniqueParams[[j]]

      if (length(dim(param_array)) == 3 && dim(param_array)[3] >= state_idx) {
        # Extract slice for this state
        param_slice <- param_array[, , state_idx, drop = FALSE]
        params[[i]][[j]] <- param_slice
      } else if (is.numeric(param_array) && length(param_array) >= state_idx) {
        # Handle vector parameters
        params[[i]][[j]] <- array(param_array[state_idx], dim = c(1, 1, 1))
      } else {
        # Default value
        params[[i]][[j]] <- array(0, dim = c(1, 1, 1))
      }
    }

    # For normal distribution, ensure proper naming
    if (inherits(mdobj, "normal")) {
      names(params[[i]]) <- c("mu", "sigma")
    }
  }

  dp <- list()

  dp$data <- x
  dp$n <- nrow(x)
  dp$mixingDistribution <- mdobj
  dp$states <- states
  dp$uniqueParams <- uniqueParams
  dp$params <- params
  dp$alpha <- alpha
  dp$beta <- beta

  class(dp) <- append(class(dp), c("markov", "dirichletprocess", class(mdobj)[-1]))

  # Set cpp preference for this object
  if (cpp) {
    options(dirichletprocess.use_cpp = TRUE)
  } else {
    options(dirichletprocess.use_cpp = FALSE)
  }

  return(dp)
}
