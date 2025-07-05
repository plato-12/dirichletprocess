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
  parameter_samples <- list(mu = numeric(no_draws), nu = numeric(no_draws))
  accept_count <- 0

  # Convert data to matrix if needed
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  }

  # Extract current values
  current_mu <- as.numeric(current_params$mu)
  current_nu <- as.numeric(current_params$nu)

  # Store first sample
  parameter_samples$mu[1] <- current_mu
  parameter_samples$nu[1] <- current_nu

  # Current log-likelihood and prior
  current_theta <- list(mu = array(current_mu, c(1,1,1)),
                        nu = array(current_nu, c(1,1,1)))
  current_log_lik <- sum(log(pmax(Likelihood(mixingDistribution, x[,1], current_theta), 1e-300)))
  current_log_prior <- log(max(PriorDensity(mixingDistribution, current_theta), 1e-300))

  # MCMC loop
  for (i in 2:no_draws) {
    # Propose new parameters
    proposed_params <- MhParameterProposal(mixingDistribution, current_theta)
    proposed_mu <- as.numeric(proposed_params$mu)
    proposed_nu <- as.numeric(proposed_params$nu)

    # Calculate proposed log-likelihood and prior
    proposed_log_lik <- sum(log(pmax(Likelihood(mixingDistribution, x[,1], proposed_params), 1e-300)))
    proposed_log_prior <- log(max(PriorDensity(mixingDistribution, proposed_params), 1e-300))

    # Calculate acceptance ratio
    log_ratio <- (proposed_log_lik + proposed_log_prior) -
      (current_log_lik + current_log_prior)

    # Accept or reject
    if (log(runif(1)) < log_ratio) {
      current_mu <- proposed_mu
      current_nu <- proposed_nu
      current_theta <- proposed_params
      current_log_lik <- proposed_log_lik
      current_log_prior <- proposed_log_prior
      accept_count <- accept_count + 1
    }

    # Store sample
    parameter_samples$mu[i] <- current_mu
    parameter_samples$nu[i] <- current_nu
  }

  return(list(
    parameter_samples = parameter_samples,
    accept_ratio = accept_count / (no_draws - 1)
  ))
}

#' @export
MetropolisHastings.default <- MetropolisHastings.beta
