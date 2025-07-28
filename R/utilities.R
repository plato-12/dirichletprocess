#' Generate a weighted function.
#'
#' @param func Function that is used of the form func(x, params).
#' @param weights Weighting of each cluster.
#' @param params Cluster parameter list
#' @return weighted function
#'
#' @export
weighted_function_generator <- function(func, weights, params) {

  weights <- weights/sum(weights)

  weightedFunc <- function(y) {

    if (is.matrix(y) || is.data.frame(y)){
      out <- numeric(nrow(y))
      y <- as.matrix(y)
    }
    else{
      out <- numeric(length(y))
    }
    cumWeight <- 0
    for (i in seq_along(weights)) {
      if (cumWeight > (1 - 1e-6)){
        break
      }
      
      # Handle MVNormal2 named parameter structure vs indexed structure
      if (!is.null(names(params)) && all(c("mu", "sig") %in% names(params))) {
        # MVNormal2 case: named parameters (mu, sig)
        cl_params <- list(
          mu = params$mu[, , i, drop = FALSE],
          sig = params$sig[, , i, drop = FALSE]
        )
      } else {
        # Standard case: indexed parameters - preserve names from original structure
        cl_params <- vector("list", length = length(params))
        param_names <- names(params)
        
        for (j in seq_along(params)) {
          param_val <- params[[j]][, , i, drop = FALSE]
          # Convert array format to scalar/vector for likelihood functions
          if (length(dim(param_val)) > 0) {
            # Remove array dimensions for scalar parameters
            if (all(dim(param_val) == c(1, 1, 1))) {
              cl_params[[j]] <- as.numeric(param_val)
            } else if (prod(dim(param_val)) == length(param_val)) {
              # For other cases, convert to vector
              cl_params[[j]] <- as.numeric(param_val)
            } else {
              cl_params[[j]] <- param_val
            }
          } else {
            cl_params[[j]] <- param_val
          }
        }
        
        # Preserve parameter names from original structure
        if (!is.null(param_names)) {
          names(cl_params) <- param_names
        }
        
        # Validate that we have sufficient parameters for likelihood functions
        if (length(cl_params) < 2) {
          # Skip this iteration if insufficient parameters
          next
        }
      }
      
      # Additional validation: ensure parameters are not empty or invalid
      if (any(sapply(cl_params, function(p) is.null(p) || length(p) == 0))) {
        next
      }
      
      out <- out + weights[i] * func(y, cl_params)
      cumWeight <- cumWeight + weights[i]
    }

    return(out)
  }
  return(weightedFunc)
}

dpareto <- function(x, xm, alpha) ifelse(x > xm, alpha * xm^alpha/(x^(alpha + 1)),
  0)
ppareto <- function(q, xm, alpha) ifelse(q > xm, 1 - (xm/q)^alpha, 0)
qpareto <- function(p, xm, alpha) ifelse(p < 0 | p > 1, NaN, xm * (1 - p)^(-1/alpha))
rpareto <- function(n, xm, alpha) qpareto(runif(n), xm, alpha)

VectorToArray <- function(paramVector){

  paramList <- vector("list", length(paramVector))

  paramList <- lapply(seq_along(paramList), function(i) array(paramVector[i], dim = c(1,1,1)))

  return(paramList)
}




