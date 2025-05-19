# Metropolis-Hastings Components Analysis

## Algorithm Overview

These components form the core of posterior sampling for non-conjugate Dirichlet Process Mixture Models. They are used within `PosteriorDraw()` for non-conjugate kernels to sample from the posterior distribution of cluster parameters $p(\theta | x, G_0)$, where `x` is the data assigned to a particular cluster.

The main functions involved are:
- `MetropolisHastings()` (S3 generic with methods like `MetropolisHastings.default` and `MetropolisHastings.weibull`)
- `MhParameterProposal()` (S3 generic, e.g., `MhParameterProposal.beta`, `MhParameterProposal.weibull`)
- `PriorDensity()` (S3 generic, e.g., `PriorDensity.beta`, `PriorDensity.weibull`)
- `Likelihood()` (S3 generic, already analyzed separately, but crucial here)

## `MetropolisHastings.*` Methods

### Algorithm Overview (`MetropolisHastings.default` from `metropolis_hastings.R`)
1. Initialize: Set `parameter_samples` array. The first sample is the `start_pos`.
2. Calculate initial log-posterior components: `old_prior = log(PriorDensity(mixingDistribution, old_param))` and `old_Likelihood = sum(log(Likelihood(mixingDistribution, x, old_param)))`.
3. For `i` from 1 to `no_draws - 1`:
   a. Propose new parameters: `prop_param = MhParameterProposal(mixingDistribution, old_param)`.
   b. Calculate log-posterior components for the proposal: `new_prior = log(PriorDensity(mixingDistribution, prop_param))` and `new_Likelihood = sum(log(Likelihood(mixingDistribution, x, prop_param)))`.
   c. Calculate acceptance probability `accept_prob = min(1, exp(new_prior + new_Likelihood - old_prior - old_Likelihood))`.
      - Handle `NA` or invalid `accept_prob` by setting it to 0.
   d. If `runif(1) < accept_prob` (accept proposal):
      - `sampled_param = prop_param`
      - `old_Likelihood = new_Likelihood`
      - `old_prior = new_prior`
   e. Else (reject proposal):
      - `sampled_param = old_param`
   f. Store `sampled_param` in `parameter_samples`.
   g. Update `old_param = sampled_param`.
4. Return `list(parameter_samples, accept_ratio)`.

### `MetropolisHastings.weibull` (metropolis_hastings.R)
- This is a specialized version that incorporates a Gibbs step for one of the parameters.
- After proposing the shape parameter `alpha` (from `MhParameterProposal.weibull`), it directly samples the scale parameter `lambda` from its full conditional posterior `Gamma(length(x) + priorParameters[2], sum(x^alpha_prop) + priorParameters[3])`, (actually sampling $1/\lambda$ then inverting to get $\lambda$).
- The acceptance probability is then based on the joint proposal (MH for `alpha`, Gibbs for `lambda`).

## `MhParameterProposal.*` Methods

### Algorithm Overview
These S3 methods define how a new candidate parameter set `new_params` is proposed based on the `old_params` and `mdObj$mhStepSize`.
- Typically involve adding random noise (e.g., from a Normal distribution scaled by `mhStepSize`) to the `old_params`.
- Must respect parameter constraints (e.g., positivity for scale/variance parameters, bounds for probabilities).

### Specific Implementations
- `MhParameterProposal.beta` (beta_uniform_gamma.R): Adds Normal noise to `mu` (mean), clips to `[0, maxT]`. Adds Normal noise to `nu` (scale/precision) and takes absolute value to ensure positivity.
- `MhParameterProposal.weibull` (weibull_uniform_gamma.R): Adds Normal noise to shape `alpha` and takes absolute value.
- `MhParameterProposal.beta2` (beta_uniform_pareto.R): Similar to `.beta` but for potentially different step sizes or proposal logic if needed.

## `PriorDensity.*` Methods

### Algorithm Overview
These S3 methods calculate the log-prior density $\log p(\theta)$ for a given parameter set `theta` based on the `mdObj$priorParameters`.

### Specific Implementations
- `PriorDensity.beta` (beta_uniform_gamma.R): `dunif(mu, 0, maxT) * dgamma(1/nu, shape, rate)` (product of densities for mu and nu).
- `PriorDensity.weibull` (weibull_uniform_gamma.R): `dunif(shape_a, 0, phi) * dgamma(1/scale_b, alpha0, beta0)`. (The R code only shows `dunif(theta[[1]], 0, priorParameters[1])` suggesting the scale part might be integrated out or handled differently in the Weibull MH or it's a simplified prior density for MH purposes, which is common if one parameter has a direct sampler).
- `PriorDensity.beta2` (beta_uniform_pareto.R): `dunif(mu, 0, maxT) * dpareto(nu, muLim, gamma_shape)`.

## Key Data Structures Used

- `mixingDistribution` (`mdObj`): Contains prior parameters, `mhStepSize`, and kernel type.
- `x`: Data points for likelihood calculation.
- `start_pos`: Initial parameter values for the MH chain.
- `no_draws`: Number of MH iterations.
- `old_param`, `prop_param`, `sampled_param`: Parameter lists/vectors.

## Performance Considerations

1.  **Number of Draws (`no_draws`):** The main loop runs this many times. This is often set by `dpObj$mhDraws`.
2.  **`Likelihood()` Call:** This is the most expensive part within each MH iteration, especially if the number of data points `x` in the cluster is large.
3.  **`PriorDensity()` Call:** Usually less expensive than `Likelihood()`.
4.  **`MhParameterProposal()` Call:** Typically very fast (random draws and arithmetic).
5.  **Log/Exp and Arithmetic:** Standard arithmetic operations.
6.  **Acceptance Rate:** Low acceptance rates mean many proposals (and their likelihood/prior calculations) are wasted, leading to inefficiency and poor mixing.

## Dependencies

- `Likelihood.*()` methods for the specific kernel.
- `PriorDensity.*()` methods for the specific kernel.
- `MhParameterProposal.*()` methods for the specific kernel.
- Base R random number generators (e.g., `rnorm`, `runif`).

