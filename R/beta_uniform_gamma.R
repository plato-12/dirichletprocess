#' Create a Beta mixing distribution.
#'
#' See \code{\link{DirichletProcessBeta}} for the default prior and hyper prior distributions.
#'
#' @param priorParameters The prior parameters for the base measure.
#' @param mhStepSize The Metropolis Hastings step size. A numeric vector of length 2.
#' @param maxT The upper bound of the Beta distribution. Defaults to 1 for the standard Beta distribution.
#' @param hyperPriorParameters The parameters for the hyper prior.
#' @return A mixing distribution object.
#' @export
BetaMixtureCreate <- function(priorParameters = c(2, 8), mhStepSize = c(1, 1), maxT = 1,
                              hyperPriorParameters = c(1, 0.125)) {

  mdObj <- MixingDistribution("beta",
                              priorParameters, "nonconjugate",
                              mhStepSize, hyperPriorParameters)
  mdObj$maxT <- maxT
  return(mdObj)
}

#' @export
#' @rdname Likelihood
Likelihood.beta <- function(mdObj, x, theta) {
  maxT <- mdObj$maxT
  x <- as.numeric(x)

  # Validate theta structure
  if (!is.list(theta)) {
    stop("theta must be a list with mu and nu components")
  }

  if (!all(c("mu", "nu") %in% names(theta))) {
    stop("theta must contain 'mu' and 'nu' components")
  }

  # Extract parameters with proper handling for various formats
  mu <- if (is.array(theta$mu)) {
    # Handle 3D arrays (dim = c(1,1,n))
    if (length(dim(theta$mu)) == 3) {
      as.numeric(theta$mu[,,, drop = TRUE])
    } else {
      as.numeric(theta$mu)
    }
  } else if (is.list(theta$mu)) {
    unlist(theta$mu)
  } else {
    as.numeric(theta$mu)
  }

  nu <- if (is.array(theta$nu)) {
    # Handle 3D arrays (dim = c(1,1,n))
    if (length(dim(theta$nu)) == 3) {
      as.numeric(theta$nu[,,, drop = TRUE])
    } else {
      as.numeric(theta$nu)
    }
  } else if (is.list(theta$nu)) {
    unlist(theta$nu)
  } else {
    as.numeric(theta$nu)
  }

  # Ensure we have valid values
  mu <- mu[!is.na(mu)]
  nu <- nu[!is.na(nu)]

  if (length(mu) == 0 || length(nu) == 0) {
    return(rep(1e-300, length(x)))
  }

  # Ensure mu and nu have the same length
  n_params <- max(length(mu), length(nu))
  if (length(mu) == 1 && n_params > 1) {
    mu <- rep(mu, n_params)
  }
  if (length(nu) == 1 && n_params > 1) {
    nu <- rep(nu, n_params)
  }

  # Calculate likelihood
  n_clusters <- length(mu)
  if (length(x) == 1) {
    # Single observation
    lik <- numeric(n_clusters)
    for (k in 1:n_clusters) {
      if (mu[k] > 0 && mu[k] < maxT && nu[k] > 0) {
        a <- (mu[k] * nu[k]) / maxT
        b <- (1 - mu[k]/maxT) * nu[k]

        if (a > 0 && b > 0 && x >= 0 && x <= maxT) {
          lik[k] <- (1/maxT) * dbeta(x/maxT, a, b)
        } else {
          lik[k] <- 1e-300
        }
      } else {
        lik[k] <- 1e-300
      }
    }
    return(if (n_clusters == 1) lik[1] else lik)
  } else {
    # Multiple observations - return matrix
    lik <- matrix(1e-300, nrow = length(x), ncol = n_clusters)
    for (k in 1:n_clusters) {
      if (mu[k] > 0 && mu[k] < maxT && nu[k] > 0) {
        a <- (mu[k] * nu[k]) / maxT
        b <- (1 - mu[k]/maxT) * nu[k]

        if (a > 0 && b > 0) {
          valid_idx <- x >= 0 & x <= maxT
          lik[valid_idx, k] <- (1/maxT) * dbeta(x[valid_idx]/maxT, a, b)
        }
      }
    }
    return(lik)
  }
}

#' @export
#' @rdname PriorDraw
PriorDraw.beta <- function(mdObj, n = 1) {

  priorParameters <- mdObj$priorParameters
  mu <- runif(n, 0, mdObj$maxT)
  
  # Draw gamma values and handle potential NAs
  gamma_values <- rgamma(n, shape = priorParameters[1], rate = priorParameters[2])
  
  # Handle NA values that can occur with extreme parameters
  if (any(is.na(gamma_values))) {
    gamma_values[is.na(gamma_values)] <- 1.0  # Default to reasonable value
  }
  
  # Ensure we don't divide by zero
  gamma_values[gamma_values == 0] <- 1e-04
  
  nu <- 1/gamma_values

  theta <- list(mu = array(mu, c(1, 1, n)), nu = array(nu, c(1, 1, n)))
  return(theta)
}

#' @export
#' @rdname PriorDensity
PriorDensity.beta <- function(mdObj, theta) {

  priorParameters <- mdObj$priorParameters
  mu <- theta[[1]]
  nu <- theta[[2]]

  muDensity <- dunif(mu, 0, mdObj$maxT)

  nuDensity <- dgamma(1/nu, priorParameters[1], priorParameters[2]) * (1/nu^2)

  if(is.infinite(nuDensity) | is.na(nuDensity)){
    nuDensity <- 1e-10 # Return a very small number instead of Inf or NA
  }

  thetaDensity <- muDensity * nuDensity
  return(as.numeric(thetaDensity))
}

#' @export
#' @rdname PriorParametersUpdate
PriorParametersUpdate.beta <- function(mdObj, clusterParameters, n = 1) {

  hyperPriorParameters <- mdObj$hyperPriorParameters
  priorParameters <- mdObj$priorParameters

  numClusters <- dim(clusterParameters[[1]])[3]

  posteriorShape <- hyperPriorParameters[1] + priorParameters[1] * numClusters
  posteriorRate <- hyperPriorParameters[2] + sum(1/clusterParameters[[2]])

  newGamma <- rgamma(n, posteriorShape, posteriorRate)

  newPriorParameters <- matrix(c(priorParameters[1], newGamma), nrow = 1, ncol = 2)
  mdObj$priorParameters <- newPriorParameters

  return(mdObj)
}

#' @export
#' @rdname MhParameterProposal
MhParameterProposal.beta <- function(mdObj, old_params) {

  mhStepSize <- mdObj$mhStepSize

  new_params <- old_params

  # Extract current values
  old_mu <- as.numeric(old_params[[1]])
  old_nu <- as.numeric(old_params[[2]])

  # Propose new mu
  new_mu <- old_mu + mhStepSize[1] * rnorm(1, 0, 2.4)
  if (new_mu > mdObj$maxT || new_mu < 0) {
    new_mu <- old_mu
  }

  # Propose new nu (ensure positive)
  new_nu <- abs(old_nu + mhStepSize[2] * rnorm(1, 0, 2.4))
  
  # Handle NA values and ensure minimum values
  if (is.na(new_nu) || new_nu == 0) {
    new_nu <- 1e-04
  }
  
  if (is.na(new_mu)) {
    new_mu <- old_mu
  }

  # Return in proper format
  new_params[[1]] <- array(new_mu, dim = c(1, 1, 1))
  new_params[[2]] <- array(new_nu, dim = c(1, 1, 1))

  return(new_params)
}

#' @export
#' @rdname PenalisedLikelihood
PenalisedLikelihood.beta <- function(mdObj, x) {
  if (length(x) == 0) {
    return(PriorDraw(mdObj, 1))
  }

  x <- as.numeric(x)
  x <- x[x > 0 & x < mdObj$maxT]  # Remove boundary values

  if (length(x) == 0) {
    return(PriorDraw(mdObj, 1))
  }

  # Method of moments estimation
  x_norm <- x / mdObj$maxT
  x_mean <- mean(x_norm)
  x_var <- var(x_norm)

  # Handle edge cases
  if (is.na(x_var) || x_var < 1e-10) {
    x_var <- 0.01
  }

  if (x_mean <= 0.01) x_mean <- 0.01
  if (x_mean >= 0.99) x_mean <- 0.99

  # Calculate parameters
  common <- x_mean * (1 - x_mean) / x_var - 1
  if (common <= 0) {
    # Fallback to prior
    return(PriorDraw(mdObj, 1))
  }

  mu_est <- x_mean * mdObj$maxT
  tau_est <- common

  # Return in the expected format
  return(list(
    mu = array(mu_est, dim = c(1, 1, 1)),
    nu = array(tau_est, dim = c(1, 1, 1))
  ))
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
