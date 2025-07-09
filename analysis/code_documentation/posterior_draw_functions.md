# PosteriorDraw Functions (`PosteriorDraw.*`) Analysis

## Algorithm Overview

The `PosteriorDraw()` S3 generic methods are responsible for drawing `n` samples of parameters `theta` from the posterior distribution $p(\theta | x, G_0)$, given the data `x` assigned to a specific cluster and the base measure `G_0` (via `mdObj`). This is central to the `ClusterParameterUpdate` step.

## Implementation Structure

S3 dispatch is used based on the class of `mdObj`, distinguishing mainly between conjugate and non-conjugate mixing distributions:
- `PosteriorDraw.conjugate` (implicit, specific methods for each conjugate kernel):
  - `PosteriorDraw.normal()` (normal_inverse_gamma.R): Samples from Normal-Inverse Gamma posterior.
  - `PosteriorDraw.normalFixedVariance()` (normal_fixed_variance.R): Samples from Normal posterior for the mean.
  - `PosteriorDraw.exponential()` (exponential_gamma.R): Samples from Gamma posterior for the rate.
  - `PosteriorDraw.mvnormal()` (mvnormal_normal_wishart.R): Samples from Normal-Wishart posterior.
- `PosteriorDraw.nonconjugate` (mixing_distribution_posterior_draw.R): Generic for non-conjugate. This method typically wraps `MetropolisHastings()`.
  - Specific non-conjugate kernels like `PosteriorDraw.beta()` (beta_uniform_gamma.R), `PosteriorDraw.weibull()` (weibull_uniform_gamma.R), `PosteriorDraw.mvnormal2()` (mvnormal_semi_conjugate.R) implement this by calling `MetropolisHastings` or having their own MH-like loop.

## Algorithmic Steps (Conjugate Cases)

1. Calculate posterior hyperparameters: The function first calls `PosteriorParameters(mdObj, x)` to compute the parameters of the analytical posterior distribution based on the prior parameters from `mdObj` and the data `x` in the current cluster.
2. Draw `n` samples from this analytical posterior distribution using base R random number generators.
   - E.g., for `PosteriorDraw.normal()`: Draw `lambda` (precision) from `Gamma(alpha_n, beta_n)`, then `mu` from `Normal(mu_n, (kappa_n * lambda)^-0.5)`.
3. Structure the samples into the standard list of arrays format.

## Algorithmic Steps (Non-conjugate Cases - General via `PosteriorDraw.nonconjugate`)

1. Determine starting position for Metropolis-Hastings: 
   - If `start_pos` is provided via `...`, use it.
   - Else, call `PenalisedLikelihood(mdObj, x)` to find a good starting point (e.g., mode of a penalized likelihood) or fall back to `PriorDraw(mdObj, 1)` if `PenalisedLikelihood` is not implemented for the specific `mdObj`.
2. Call `MetropolisHastings(mdObj, x, start_pos, no_draws = n)` to obtain `n` samples from the posterior via MCMC.
3. Format the output samples from `MetropolisHastings` into the standard list of arrays.

### Specific Non-conjugate Implementations (e.g., `PosteriorDraw.weibull`):
- May have custom logic within the MH loop, for example, if one parameter can be updated via Gibbs sampling conditional on the other MH-sampled parameter (as in Weibull where `lambda` is sampled from its conditional posterior given `alpha`).
- `PosteriorDraw.weibull` explicitly samples `lambda` from its conditional posterior (an Inverse Gamma derived from the Gamma on $1/\lambda$) given the current `alpha` proposed by `MhParameterProposal.weibull`.

## Key Data Structures Used

- `mdObj`: Mixing distribution object, containing prior parameters and methods like `Likelihood`, `PriorDensity`, `MhParameterProposal`.
- `x`: Data points assigned to the current cluster.
- `n`: Number of posterior samples to draw.
- `...` (often `start_pos`): Starting parameters for non-conjugate MCMC.
- Output: A list of arrays, similar to `PriorDraw` output.

## Performance Considerations

1.  **Conjugate Models**: Performance is dictated by the cost of calculating posterior hyperparameters (usually involves summing over data `x`) and then drawing `n` samples. Generally efficient.
2.  **Non-conjugate Models**: Significantly more computationally intensive.
    -   Dominated by the `MetropolisHastings` call, which runs for `n` (often `mhDraws`) iterations.
    -   Each MH iteration involves: parameter proposal, `Likelihood()` call (summed over data `x`), `PriorDensity()` call, and acceptance step.
    -   The cost of `Likelihood(mdObj, x, theta)` is crucial here as `x` can be many data points.
    -   `PenalisedLikelihood()` (if used for `start_pos`) can involve an optimization routine which adds to the setup cost.

## Dependencies

- For conjugate cases: `PosteriorParameters()` method, base R RNGs.
- For non-conjugate cases: `MetropolisHastings()` generic (and its specific methods), `MhParameterProposal()`, `Likelihood()`, `PriorDensity()`, and potentially `PenalisedLikelihood()` or `PriorDraw()` for starting positions.

