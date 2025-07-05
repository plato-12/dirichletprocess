MetropolisHastings <- function(mixingDistribution, x, start_pos, no_draws=100){
  UseMethod("MetropolisHastings", mixingDistribution)
}

MetropolisHastings.default <- function(mixingDistribution, x, start_pos, no_draws = 100) {
  parameter_samples <- vector("list", length(start_pos))
  for (i in seq_along(start_pos)) {
    parameter_samples[[i]] <- array(dim = c(dim(start_pos[[i]])[1:2], no_draws))
    parameter_samples[[i]][, , 1] <- start_pos[[i]][, , 1]
  }

  accept_count <- 0
  old_param <- start_pos

  old_prior <- log(PriorDensity(mixingDistribution, old_param))
  old_Likelihood <- sum(log(Likelihood(mixingDistribution, x, old_param)))

  for (i in seq_len(no_draws - 1)) {

    prop_param <- MhParameterProposal(mixingDistribution, old_param)

    new_prior <- log(PriorDensity(mixingDistribution, prop_param))
    new_Likelihood <- sum(log(Likelihood(mixingDistribution, x, prop_param)))

    accept_prob <- min(1, exp(new_prior + new_Likelihood - old_prior - old_Likelihood))

    if (is.na(accept_prob) || !is.finite(accept_prob) ) {
      accept_prob <- 0
    }

    if (runif(1) < accept_prob) {
      accept_count <- accept_count + 1
      sampled_param <- prop_param
      old_Likelihood <- new_Likelihood
      old_prior <- new_prior
    } else {
      sampled_param <- old_param
    }

    old_param <- sampled_param
    for (j in seq_along(start_pos)) {
      parameter_samples[[j]][, , i + 1] <- sampled_param[[j]]
    }
  }

  accept_ratio <- accept_count/no_draws

  return(list(parameter_samples = parameter_samples, accept_ratio = accept_ratio))
}

MetropolisHastings.weibull <- function(mixingDistribution, x, start_pos, no_draws=100){

  lamSamp <- 1/rgamma(1, length(x)+mixingDistribution$priorParameters[2],
                      sum(x^c(start_pos[[1]])) + mixingDistribution$priorParameters[3])
  start_pos[[2]] <- array(lamSamp, dim=c(1,1,1))

  parameter_samples <- list()
  for (i in seq_along(start_pos)) {
    parameter_samples[[i]] <- array(dim = c(dim(start_pos[[i]])[1:2], no_draws))
    parameter_samples[[i]][, , 1] <- start_pos[[i]][, , 1]
  }

  accept_count <- 0
  old_param <- start_pos

  old_prior <- log(PriorDensity(mixingDistribution, old_param))
  old_Likelihood <- sum(log(Likelihood(mixingDistribution, x, old_param)))

  for(i in seq_len(no_draws-1)){

    prop_param <- MhParameterProposal(mixingDistribution, old_param)
    lamSamp <- 1/rgamma(1, length(x)+mixingDistribution$priorParameters[2],
                        sum(x^c(prop_param[[1]])) + mixingDistribution$priorParameters[3])
    prop_param[[2]] <- array(lamSamp, dim=c(1,1,1))

    new_prior <- log(PriorDensity(mixingDistribution, prop_param))
    new_Likelihood <- sum(log(Likelihood(mixingDistribution, x, prop_param)))

    accept_prob <- min(1, exp(new_prior + new_Likelihood - old_prior - old_Likelihood))

    if (is.na(accept_prob)) {
      accept_prob <- 0
    }

    if (runif(1) < accept_prob) {
      accept_count <- accept_count + 1
      sampled_param <- prop_param
      old_Likelihood <- new_Likelihood
      old_prior <- new_prior
    } else {
      sampled_param <- old_param
    }

    old_param <- sampled_param
    for (j in seq_along(start_pos)) {
      parameter_samples[[j]][, , i + 1] <- sampled_param[[j]]
    }

  }
  accept_ratio <- accept_count/no_draws

  return(list(parameter_samples = parameter_samples, accept_ratio = accept_ratio))

}

#' @export
MetropolisHastings.beta <- function(mixingDistribution, x, start_pos, no_draws) {
  # Initialize
  current_params <- start_pos
  parameter_samples <- list()

  # Initialize parameter sample arrays following the pattern from default method
  for (i in seq_along(start_pos)) {
    parameter_samples[[i]] <- array(dim = c(dim(start_pos[[i]])[1:2], no_draws))
    parameter_samples[[i]][, , 1] <- start_pos[[i]][, , 1]
  }

  accept_count <- 0

  # Convert data to matrix if needed
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }

  # Current parameters
  old_param <- start_pos

  # Current log-likelihood and prior
  old_prior <- log(PriorDensity(mixingDistribution, old_param))
  old_Likelihood <- sum(log(Likelihood(mixingDistribution, x[,1], old_param)))

  # MCMC loop
  for (i in seq_len(no_draws - 1)) {
    # Propose new parameters
    prop_param <- MhParameterProposal(mixingDistribution, old_param)

    # Calculate proposed log-likelihood and prior
    new_prior <- log(PriorDensity(mixingDistribution, prop_param))
    new_Likelihood <- sum(log(Likelihood(mixingDistribution, x[,1], prop_param)))

    # Calculate acceptance ratio
    accept_prob <- min(1, exp(new_prior + new_Likelihood - old_prior - old_Likelihood))

    if (is.na(accept_prob) || !is.finite(accept_prob)) {
      accept_prob <- 0
    }

    # Accept or reject
    if (runif(1) < accept_prob) {
      accept_count <- accept_count + 1
      sampled_param <- prop_param
      old_Likelihood <- new_Likelihood
      old_prior <- new_prior
    } else {
      sampled_param <- old_param
    }

    # Update old_param for next iteration
    old_param <- sampled_param

    # Store sample
    for (j in seq_along(start_pos)) {
      parameter_samples[[j]][, , i + 1] <- sampled_param[[j]]
    }
  }

  accept_ratio <- accept_count / no_draws

  return(list(
    parameter_samples = parameter_samples,
    accept_ratio = accept_ratio
  ))
}
