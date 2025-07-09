# PriorDraw Functions (`PriorDraw.*`) Analysis

## Algorithm Overview

The `PriorDraw()` S3 generic methods are responsible for drawing `n` samples of parameters `theta` from the base measure `G_0` associated with a specific mixing distribution `mdObj`. These draws are used to parameterize new clusters when they are formed, and to generate auxiliary parameters in Neal's Algorithm 8 for non-conjugate mixtures.

## Implementation Structure

S3 dispatch is used based on the class of `mdObj`. Each mixing distribution has its own method defining how to sample from its specific base measure `G_0`:
- `PriorDraw.normal()`: Normal-Inverse Gamma base measure for Gaussian kernel.
- `PriorDraw.normalFixedVariance()`: Normal base measure for mean (variance is fixed).
- `PriorDraw.beta()`: Uniform for mean `mu`, Inverse-Gamma for scale `nu`.
- `PriorDraw.beta2()`: Uniform for mean `mu`, Pareto for scale `nu`.
- `PriorDraw.weibull()`: Uniform for shape `a`, Inverse-Gamma for scale `lambda` (used as `b`).
- `PriorDraw.exponential()`: Gamma base measure for the rate parameter.
- `PriorDraw.mvnormal()`: Normal for mean `mu`, Wishart for precision matrix `sig` (likely covariance in practice).
- `PriorDraw.mvnormal2()`: Normal for mean `mu`, Inverse-Wishart for covariance matrix `sig`.
- `PriorDraw.hierarchical()`: Samples from the global discrete measure $G_0$ which itself is a draw from a DP. It samples an index based on global stick-breaking weights `pi_k` and returns the corresponding global parameter `theta_k`.

## Algorithmic Steps (General)

1. Retrieve prior parameters from `mdObj$priorParameters`.
2. For each of the `n` samples requested:
   a. Draw each component of the parameter vector `theta` from its specified prior distribution using base R's random number generators (e.g., `rnorm()`, `rgamma()`, `runif()`, `rWishart()`).
   b. The parameters drawn might be independent or dependent (e.g., in Normal-Inverse Gamma, the variance is drawn first, then the mean conditional on the variance).
3. Structure the `n` samples into a list of arrays. Each list element corresponds to a parameter (e.g., `theta[[1]]` for $\mu$, `theta[[2]]` for $\sigma$), and each array has dimensions appropriate for the parameter (e.g., `1x1xn` for univariate scalar parameters, `dx1xn` for mean vectors, `dxdxn` for covariance matrices).

## Specific Kernel Details

### `PriorDraw.normal` (normal_inverse_gamma.R)
- Draws `lambda` (precision) from `Gamma(alpha0, beta0)`.
- Draws `mu` (mean) from `Normal(mu0, (kappa0 * lambda)^-0.5)`.
- Returns `list(array(mu), array(sqrt(1/lambda)))`.

### `PriorDraw.beta` (beta_uniform_gamma.R)
- Draws `mu` from `Uniform(0, mdObj$maxT)`.
- Draws `nu` (scale) from `1 / Gamma(alpha0, beta0)` (i.e., Inverse-Gamma).

### `PriorDraw.weibull` (weibull_uniform_gamma.R)
- Draws shape `a` from `Uniform(0, priorParameters[1])` (phi).
- Draws scale `lambda` (representing `b`) from `1 / Gamma(priorParameters[2], priorParameters[3])` (alpha0, beta0 for Inv-Gamma).

### `PriorDraw.mvnormal` (mvnormal_normal_wishart.R)
- Draws covariance matrix `sig` from `Wishart(nu, Lambda)` (where `Lambda` is the scale matrix for `rWishart`).
- Draws mean `mu` from `Normal(mu0, solve(sig * kappa0))` (multivariate normal).

### `PriorDraw.hierarchical` (mixing_distribution_prior_draw.R)
- Uses probabilities `mdObj$pi_k` (global stick-breaking weights for $G_0$).
- Samples `n` indices `ind` based on these probabilities.
- Returns `lapply(mdObj$theta_k, function(x) x[,, ind, drop = FALSE])`, effectively selecting `n` parameters from the existing discrete global atoms `theta_k`.

## Key Data Structures Used

- `mdObj$priorParameters`: Contains the hyperparameters for the base measure `G_0`.
- `n`: The number of parameter sets to draw.
- Output: A list of arrays, where each array holds `n` samples of a specific parameter.

## Performance Considerations

1.  **Efficiency of RNGs**: Performance depends on the efficiency of R's underlying random number generators for distributions like Gamma, Normal, Uniform, Wishart.
2.  **Looping for `n` samples**: Each of the `n` samples is typically drawn independently.
3.  **Complexity of Prior**: More complex prior structures (e.g., hierarchical dependencies between parameters) will take longer to sample.
4.  **`PriorDraw.hierarchical`**: Depends on sampling from a discrete distribution, which is generally efficient.

## Dependencies

- Base R random number generation functions (`rgamma`, `rnorm`, `runif`, etc.).
- `rWishart` from the `stats` package for multivariate normal base measures.

