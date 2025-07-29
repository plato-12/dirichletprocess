#' Create a Weibull mixing distribution.
#'
#' See \code{\link{DirichletProcessWeibull}} for the default prior and hyper prior distributions.
#'
#' @param priorParameters Prior parameters for the Weibull parameters
#' @param mhStepSize Metropolis Hastings Step Size
#' @param hyperPriorParameters Parameters for the hyper-priors
#' @return A mixing distribution object.
#' @export
WeibullMixtureCreate <- function(priorParameters, mhStepSize,
                                 hyperPriorParameters = c(6, 2, 1, 0.5)) {

  mdObj <- MixingDistribution("weibull", priorParameters, "nonconjugate",
                              mhStepSize, hyperPriorParameters)
  return(mdObj)
}

#' @export
Likelihood.weibull <- function(mdObj, x, theta) {
  # as.numeric(dweibull(x, theta[[1]], theta[[2]]))
  x <- as.vector(x, "numeric")
  alpha <- theta[[1]][, , , drop = TRUE]
  lambda <- theta[[2]][, , , drop = TRUE]

  # a <- alpha
  # b <- lambda^(1/alpha)
  # b[is.infinite(b)] <- 1000000000000000
  #y <- dweibull(x, a, b, log = TRUE)

  y <- as.numeric(lambda^(-1) * alpha * x^(alpha - 1) * exp(-lambda^(-1) * x^alpha))
  y[is.infinite(lambda)] <- 0
  y[x < 0] <- 0
  return(y)
}

#' @export
#' @rdname PriorDraw
PriorDraw.weibull <- function(mdObj, n = 1, ...) {

  priorParameters <- mdObj$priorParameters

  # Draw gamma values and handle potential NAs
  gamma_values <- rgamma(n, priorParameters[2], priorParameters[3])
  
  # Handle NA values that can occur with extreme parameters
  if (any(is.na(gamma_values))) {
    gamma_values[is.na(gamma_values)] <- 1.0  # Default to reasonable value
  }
  
  # Ensure we don't divide by zero
  gamma_values[gamma_values == 0] <- 1e-04
  
  lambdas <- 1/gamma_values
  
  # Draw uniform values and handle potential NAs
  alpha_values <- runif(n, 0, priorParameters[1])
  if (any(is.na(alpha_values))) {
    alpha_values[is.na(alpha_values)] <- 1.0  # Default to reasonable value
  }
  
  theta <- list(array(alpha_values, dim = c(1, 1, n)),
                array(lambdas, dim = c(1, 1, n)))
  return(theta)
}

#' @export
#' @rdname PriorDensity
PriorDensity.weibull <- function(mdObj, theta) {

  priorParameters <- mdObj$priorParameters

  # Handle different input types (matrix or list)
  if (is.matrix(theta)) {
    theta_val <- theta[1, 1]
  } else if (is.list(theta)) {
    # Handle different parameter dimensions safely
    if (is.array(theta[[1]])) {
      param_dims <- dim(theta[[1]])
      if (length(param_dims) == 3) {
        theta_val <- theta[[1]][1,1,1]
      } else if (length(param_dims) == 2) {
        theta_val <- theta[[1]][1,1]
      } else {
        theta_val <- theta[[1]][1]
      }
    } else {
      theta_val <- theta[[1]]
    }
  } else {
    theta_val <- theta[1]
  }
  
  theta_density <- dunif(as.numeric(theta_val), 0, priorParameters[1])
  #theta_density <- thetaDensity * dgamma(1/theta[[2]], priorParameters[2], priorParameters[3])
  return(theta_density)
}

#' @export
#' @rdname PosteriorDraw
PosteriorDraw.weibull <- function(mdObj, x, n = 100, ...) {

  if (missing(...)){
    start_pos <- PriorDraw(mdObj, 1)
  } else {
    start_pos <- list(...)$start_pos
  }

  mh_result <- MetropolisHastings.weibull(mdObj, x, start_pos, no_draws = n)

  theta <- list(array(mh_result$parameter_samples[[1]], dim = c(1, 1, n)),
                array(mh_result$parameter_samples[[2]],
                      dim = c(1, 1, n)))

  return(theta)
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.weibull <- function(mdObj, clusterParameters, n = 1) {

  hyperPriorParameters <- mdObj$hyperPriorParameters
  priorParameters <- mdObj$priorParameters

  numClusters <- dim(clusterParameters[[1]])[3]

  newPhi <- rpareto(n, max(clusterParameters[[1]], hyperPriorParameters[1]),
                    hyperPriorParameters[2] + numClusters)
  newGamma <- rgamma(n, hyperPriorParameters[3] + 2 * numClusters,
                     hyperPriorParameters[4] + sum(1/clusterParameters[[2]]))

  new_priorParameters <- matrix(c(newPhi[n],
                                  priorParameters[2],
                                  newGamma[n]),
                                  ncol = 3)
  mdObj$priorParameters <- new_priorParameters
  return(mdObj)
}

#' @export
MhParameterProposal.weibull <- function(mdObj, old_params) {

  mhStepSize <- mdObj$mhStepSize
  new_params <- old_params
  
  # Extract current values
  old_alpha <- as.numeric(old_params[[1]])
  old_lambda <- as.numeric(old_params[[2]])
  
  # Propose new alpha (ensure positive)
  new_alpha <- abs(old_alpha + mhStepSize[1] * rnorm(1, 0, 1.7))
  
  # Handle NA values and ensure minimum values
  if (is.na(new_alpha) || new_alpha == 0) {
    new_alpha <- 1e-04
  }
  
  # Propose new lambda (ensure positive)
  new_lambda <- abs(old_lambda + mhStepSize[2] * rnorm(1, 0, 1.7))
  
  # Handle NA values and ensure minimum values
  if (is.na(new_lambda) || new_lambda == 0) {
    new_lambda <- 1e-04
  }
  
  # Return in proper format
  new_params[[1]] <- array(new_alpha, dim = c(1, 1, 1))
  new_params[[2]] <- array(new_lambda, dim = c(1, 1, 1))
  
  return(new_params)
}

#' @export
MetropolisHastings.beta <- function(mixingDistribution, x, start_pos, no_draws) {
  # Initialize parameter storage
  parameter_samples <- list()
  for (i in seq_along(start_pos)) {
    parameter_samples[[i]] <- array(dim = c(dim(start_pos[[i]])[1:2], no_draws))
    parameter_samples[[i]][, , 1] <- start_pos[[i]][, , 1]
  }

  accept_count <- 0
  old_param <- start_pos

  # Calculate initial log prior and likelihood
  old_prior <- log(PriorDensity(mixingDistribution, old_param))
  old_likelihood <- sum(log(Likelihood(mixingDistribution, x, old_param)))

  # MCMC loop
  for (i in seq_len(no_draws - 1)) {
    # Propose new parameters
    prop_param <- MhParameterProposal(mixingDistribution, old_param)

    # Calculate new log prior and likelihood
    new_prior <- log(PriorDensity(mixingDistribution, prop_param))
    new_likelihood <- sum(log(Likelihood(mixingDistribution, x, prop_param)))

    # Calculate acceptance probability
    log_ratio <- (new_prior + new_likelihood) - (old_prior + old_likelihood)
    accept_prob <- min(1, exp(log_ratio))

    # Handle numerical issues
    if (is.na(accept_prob) || !is.finite(accept_prob)) {
      accept_prob <- 0
    }

    # Accept or reject
    if (runif(1) < accept_prob) {
      accept_count <- accept_count + 1
      sampled_param <- prop_param
      old_likelihood <- new_likelihood
      old_prior <- new_prior
    } else {
      sampled_param <- old_param
    }

    # Store parameters
    old_param <- sampled_param
    for (j in seq_along(start_pos)) {
      parameter_samples[[j]][, , i + 1] <- sampled_param[[j]][, , 1]
    }
  }

  accept_ratio <- accept_count / no_draws

  return(list(parameter_samples = parameter_samples, accept_ratio = accept_ratio))
}

#' @export
MetropolisHastings.list <- function(mixingDistribution, x, start_pos, no_draws = 100) {
  # For list objects, dispatch based on the second class in the hierarchy
  if (length(class(mixingDistribution)) > 1) {
    dist_class <- class(mixingDistribution)[2]
    ns <- getNamespace("dirichletprocess")
    
    # Handle weibull
    if (dist_class == "weibull") {
      weibull_func <- get("MetropolisHastings.weibull", envir = ns)
      return(weibull_func(mixingDistribution, x, start_pos, no_draws))
    }
    
    # Handle beta
    if (dist_class == "beta") {
      # Call the beta method directly here
      parameter_samples <- list()
      for (i in seq_along(start_pos)) {
        parameter_samples[[i]] <- array(dim = c(dim(start_pos[[i]])[1:2], no_draws))
        parameter_samples[[i]][, , 1] <- start_pos[[i]][, , 1]
      }

      accept_count <- 0
      old_param <- start_pos

      # Calculate initial log prior and likelihood
      old_prior <- log(PriorDensity(mixingDistribution, old_param))
      old_likelihood <- sum(log(Likelihood(mixingDistribution, x, old_param)))

      # MCMC loop
      for (i in seq_len(no_draws - 1)) {
        # Propose new parameters
        prop_param <- MhParameterProposal(mixingDistribution, old_param)

        # Calculate new log prior and likelihood
        new_prior <- log(PriorDensity(mixingDistribution, prop_param))
        new_likelihood <- sum(log(Likelihood(mixingDistribution, x, prop_param)))

        # Calculate acceptance probability
        log_ratio <- (new_prior + new_likelihood) - (old_prior + old_likelihood)
        accept_prob <- min(1, exp(log_ratio))

        # Handle numerical issues
        if (is.na(accept_prob) || !is.finite(accept_prob)) {
          accept_prob <- 0
        }

        # Accept or reject
        if (runif(1) < accept_prob) {
          accept_count <- accept_count + 1
          sampled_param <- prop_param
          old_likelihood <- new_likelihood
          old_prior <- new_prior
        } else {
          sampled_param <- old_param
        }

        # Store parameters
        old_param <- sampled_param
        for (j in seq_along(start_pos)) {
          parameter_samples[[j]][, , i + 1] <- sampled_param[[j]][, , 1]
        }
      }

      accept_ratio <- accept_count / no_draws
      return(list(parameter_samples = parameter_samples, accept_ratio = accept_ratio))
    }
  }
  
  # Fall back to default method
  return(MetropolisHastings.default(mixingDistribution, x, start_pos, no_draws))
}
