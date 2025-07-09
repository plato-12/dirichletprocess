# Predictive Likelihood Functions (`Predictive.*`) Analysis

## Algorithm Overview

The `Predictive()` S3 generic methods calculate the marginal likelihood of new data point(s) `x` given the base measure `G_0` (defined by `mdObj`). This is $p(x | G_0) = \int k(x | \theta) dG_0(\theta)$. These functions are primarily used in the `ClusterComponentUpdate.conjugate()` method to determine the probability of assigning a data point to a new cluster.

## Implementation Structure

S3 dispatch is used based on the class of the mixing distribution object (`mdObj`). These are typically implemented only for conjugate mixture models, as the integral is often analytically tractable in such cases:
- `Predictive.normal()` (normal_inverse_gamma.R): For the Normal-Inverse Gamma base measure.
- `Predictive.normalFixedVariance()` (normal_fixed_variance.R): For the Normal base measure (fixed variance).
- `Predictive.exponential()` (exponential_gamma.R): For the Gamma base measure.
- `Predictive.mvnormal()` (mvnormal_normal_wishart.R): For the Normal-Wishart base measure.
- Non-conjugate kernels generally do not have a simple `Predictive()` method because the integral $\int k(x | \theta) dG_0(\theta)$ is not analytically tractable. In such cases, Neal's Algorithm 8 (using auxiliary parameters) is used instead of relying on this predictive likelihood for new cluster formation.

## Algorithmic Steps (General for Conjugate Kernels)

1. Retrieve prior parameters from `mdObj$priorParameters`.
2. For each data point `x[i]` (if `x` is a vector/matrix of multiple points):
   a. Analytically derive or use known results for the marginal likelihood (predictive distribution) of `x[i]` under the specified `mdObj` (kernel `k` and base measure `G_0`).
   b. This often involves updating the prior parameters to 'posterior' parameters as if `x[i]` was observed, and then using ratios of normalizing constants or properties of standard distributions.
      - For example, in `Predictive.normal()`, it calculates the posterior parameters as if `x[i]` was the only data point, then uses a formula involving Gamma functions and ratios of prior/posterior parameters (related to the Student-t predictive distribution).
      - For `Predictive.exponential()`, it computes $p(x_i) = \frac{\Gamma(\alpha_0+1)}{\Gamma(\alpha_0)} \frac{\beta_0^{\alpha_0}}{(\beta_0+x_i)^{\alpha_0+1}}$ (if $x_i$ is a single observation and length(x[i])=1, and sum(x[i])=x[i]). The code `alphaPost <- priorParameters[1] + length(x[i]); betaPost <- priorParameters[2] + sum(x[i])` suggests it's calculating for a set of x[i] as if they formed a single new observation, which is correct for the predictive likelihood of a new data point being integrated over the prior.
3. Store the calculated predictive likelihood for each `x[i]` in `predictiveArray`.
4. Return `predictiveArray`.

## Specific Kernel Details

### `Predictive.normal` (normal_inverse_gamma.R)
- For each data point `x[i]`, it first computes what the posterior parameters (`mu_n`, `kappa_n`, `alpha_n`, `beta_n`) would be if `x[i]` were observed (by calling `PosteriorParameters(mdObj, x[i])`).
- Then, it uses the formula: `(gamma(alpha_n)/gamma(alpha0)) * ((beta0^alpha0)/(beta_n^alpha_n)) * sqrt(kappa0/kappa_n)`.
  This corresponds to the density of a non-standardized Student-t distribution, which is the predictive distribution for a Normal likelihood with Normal-Inverse Gamma prior.

### `Predictive.exponential` (exponential_gamma.R)
- For each `x[i]`, it calculates `alphaPost = alpha0 + 1` (assuming length of x[i] is 1) and `betaPost = beta0 + x[i]`.
- The predictive likelihood is `(gamma(alphaPost)/gamma(alpha0)) * ((beta0^alpha0)/(betaPost^alphaPost))` which is the density of a Beta-prime (or Gamma-Gamma) distribution scaled appropriately.

### `Predictive.mvnormal` (mvnormal_normal_wishart.R)
- Calculates posterior parameters `post_params` for each `x[i,]`.
- The formula involves `pi^(-d/2)`, ratios of `kappa0/kappa_n`, determinants of `Lambda` (prior scale matrix) and `t_n` (posterior scale matrix), and a product of Gamma function ratios. This corresponds to the density of a multivariate Student-t distribution.

## Key Data Structures Used

- `mdObj`: The mixing distribution object, containing `priorParameters` and defining the kernel and base measure.
- `x`: The data point(s) for which the predictive likelihood is being calculated (numeric vector or matrix).
- Output: A numeric vector `predictiveArray` of the same length as `x`, containing the predictive likelihoods.

## Performance Considerations

1.  **Loop over Data Points**: The functions typically loop through each data point in `x` if multiple are provided.
2.  **Complexity of Analytical Formula**: The cost per data point depends on the complexity of the analytical formula for the predictive likelihood. This can involve Gamma functions, determinants (for MVN), square roots, and other arithmetic operations.
3.  **`PosteriorParameters()` Call**: Some implementations (like `Predictive.normal`) call `PosteriorParameters()` internally for each data point, which itself involves calculations based on that point.

## Dependencies

- Base R math functions (`gamma`, `sqrt`, `log`, `det`, `pi`).
- `PosteriorParameters()` method for the specific conjugate kernel (used by some implementations).

