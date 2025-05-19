# Likelihood Functions (`Likelihood.*`) Analysis

## Algorithm Overview

The `Likelihood()` functions are S3 generic methods responsible for calculating the likelihood of one or more data points `x` given a specific set of parameters `theta` for a particular mixing distribution (kernel) defined by `mdObj`. These are fundamental for the DPMM fitting process, used in assigning points to clusters and in updating cluster parameters (especially in non-conjugate cases).

## Implementation Structure

S3 dispatch is used based on the class of the mixing distribution object (`mdObj`). Specific methods exist for each supported kernel:
- `Likelihood.normal()`: For Gaussian mixture kernel (unknown mean and variance).
- `Likelihood.normalFixedVariance()`: For Gaussian mixture kernel with fixed variance.
- `Likelihood.beta()`: For Beta mixture kernel (parameterized by mean and precision, on `[0, maxT]`).
- `Likelihood.beta2()`: For Beta mixture kernel (alternative parameterization/prior, also on `[0, maxT]`).
- `Likelihood.weibull()`: For Weibull mixture kernel.
- `Likelihood.exponential()`: For Exponential mixture kernel.
- `Likelihood.mvnormal()`: For Multivariate Normal mixture kernel (conjugate Normal-Wishart prior structure).
- `Likelihood.mvnormal2()`: For Multivariate Normal mixture kernel (semi-conjugate, independent priors on mean and covariance).
- Custom implementations for user-defined kernels would follow this pattern.

## Algorithmic Steps (General)

For a given data point(s) `x` and parameter set `theta`:
1. Extract the relevant parameters from the `theta` list (e.g., mean and standard deviation for Normal, shape and scale for Weibull, etc.). The `theta` argument is often a list where each element is an array, allowing for vectorized calculations if multiple parameter sets are provided (e.g., for different clusters or auxiliary parameters).
2. Handle the dimensionality: `theta` parameters are often structured as 3D arrays `(param_dim1, param_dim2, num_param_sets)`. The function needs to correctly extract parameters for each set if `x` is being evaluated against multiple `theta`.
3. For each data point in `x` (if `x` is a vector or matrix of multiple points) and for each parameter set in `theta`:
   a. Call the corresponding base R density function (e.g., `dnorm()`, `dbeta()`, `dweibull()`, `mvtnorm::dmvnorm()`).
   b. Adjustments for parameterization might be needed. For example, `Likelihood.beta` converts mean/precision (mu, tau) to shape1/shape2 (a,b) and scales by `1/maxT` because `dbeta` in R is for the standard [0,1] Beta.
4. Return a numeric vector or matrix of likelihood values. If multiple parameter sets in `theta` are provided, the output is typically a vector where each element corresponds to the likelihood of `x` under one parameter set (e.g., in `Likelihood.mvnormal`).

## Specific Kernel Details

### `Likelihood.normal` (normal_inverse_gamma.R)
- Parameters in `theta`: `[[1]]` is mean (`mu`), `[[2]]` is standard deviation (`sigma`).
- Uses `dnorm(x, theta[[1]], theta[[2]])`.

### `Likelihood.beta` (beta_uniform_gamma.R)
- Parameters in `theta`: `[[1]]` is mean (`mu`), `[[2]]` is precision/scale (`tau` or `nu`).
- `maxT` is retrieved from `mdObj$maxT`.
- Converts (mu, tau) to standard Beta parameters (a,b): `a = (mu * tau)/maxT`, `b = (1 - mu/maxT) * tau`.
- Calculates likelihood as `1/maxT * dbeta(x/maxT, a, b)` to account for the `[0, maxT]` support.

### `Likelihood.weibull` (weibull_uniform_gamma.R)
- Parameters in `theta`: `[[1]]` is shape (`alpha`), `[[2]]` is scale (`lambda`). The code uses $ \lambda^{-1} \alpha x^{\alpha-1} \exp(-\lambda^{-1} x^\alpha) $ which means $\lambda$ is the scale parameter `b`.
- Handles `Inf` values for `lambda` (scale) by setting likelihood to 0.
- Sets likelihood to 0 for `x < 0`.

### `Likelihood.mvnormal` & `Likelihood.mvnormal2` (mvnormal_normal_wishart.R, mvnormal_semi_conjugate.R)
- Parameters in `theta`: `[[1]]` (or `theta$mu`) is the mean vector, `[[2]]` (or `theta$sig`) is the covariance matrix.
- Uses `mvtnorm::dmvnorm(x, mean, sigma)`.
- Loops or `vapply`s through multiple parameter sets if provided in `theta`.

## Key Data Structures Used

- `mdObj`: The mixing distribution object, which might contain fixed parameters (e.g., `mdObj$sigma` in `Likelihood.normalFixedVariance`, `mdObj$maxT` in `Likelihood.beta`).
- `x`: The data point(s) (numeric vector or matrix).
- `theta`: A list of parameter arrays. Each element of the list corresponds to a parameter of the distribution (e.g., `theta[[1]]` for $\mu$, `theta[[2]]` for $\sigma^2$). The third dimension of these arrays typically allows for evaluating likelihood against multiple parameter sets.

## Performance Considerations

1.  **Vectorization**: Efficient implementations vectorize calculations over `x` if `theta` contains a single parameter set, or over `theta` if `x` is a single data point evaluated against multiple parameter sets. The use of `vapply` in some mvnormal cases suggests iteration over parameter sets.
2.  **Base R Efficiency**: Relies on the efficiency of base R's density functions (`dnorm`, `dbeta`, etc.) or `mvtnorm::dmvnorm()`.
3.  **Parameter Transformation**: For some kernels (e.g., Beta), parameters might need transformation before calling the base R density function, adding slight overhead.
4.  **Repeated Calls**: These functions are called very frequently. Even small inefficiencies can compound.

## Dependencies

- Base R density functions (e.g., `dnorm`, `dbeta`, `dweibull`, `dexp`).
- `mvtnorm::dmvnorm` for multivariate normal distributions.

