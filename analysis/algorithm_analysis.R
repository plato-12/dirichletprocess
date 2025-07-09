# Save as analysis/algorithm_analysis.R
library(dirichletprocess)

document_algorithms <- function() {
  # Create directory for documentation if it doesn't exist
  if (!dir.exists("analysis/code_documentation")) {
    dir.create("analysis/code_documentation", recursive = TRUE)
  }

  # Analyze ClusterComponentUpdate
  analyze_cluster_component_update()

  # Analyze ClusterParameterUpdate
  analyze_cluster_parameter_update()

  # Analyze UpdateAlpha
  analyze_update_alpha()

  # Analyze Likelihood functions
  analyze_likelihood_functions()

  # Analyze PriorDraw functions
  analyze_prior_draw_functions()

  # Analyze PosteriorDraw functions
  analyze_posterior_draw_functions()

  # Analyze Metropolis-Hastings related functions
  analyze_metropolis_hastings_components()

  # Analyze Predictive Likelihood functions
  analyze_predictive_likelihood_functions()

  # Analyze ClusterLabelChange function
  analyze_cluster_label_change()

  # Analyze Hierarchical Model Specific Updates
  analyze_hierarchical_updates()
}

analyze_cluster_component_update <- function() {
  sink("analysis/code_documentation/cluster_component_update.md")

  cat("# ClusterComponentUpdate Algorithm Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("This implements Neal's Algorithm 4 (for conjugate cases) and Algorithm 8 (for non-conjugate cases) for updating cluster assignments of data points in a Dirichlet Process Mixture Model (DPMM).\n\n")

  cat("## Implementation Structure\n\n")
  cat("The function is implemented with S3 dispatch based on the `dpObj` class:\n")
  cat("- `ClusterComponentUpdate()`: Generic function.\n")
  cat("- `ClusterComponentUpdate.conjugate()`: Handles DPMMs with conjugate base measures.\n")
  cat("- `ClusterComponentUpdate.nonconjugate()`: Handles DPMMs with non-conjugate base measures.\n")
  cat("- `ClusterComponentUpdate.hierarchical()`: Handles hierarchical Dirichlet Processes.\n\n")

  cat("## Algorithmic Steps (Conjugate Case - `ClusterComponentUpdate.conjugate`)\n\n")
  cat("Refers to Neal's Algorithm 4.\n")
  cat("1. For each data point `i` from 1 to `n` (total number of data points):\n")
  cat("   a. Temporarily remove data point `y[i,]` from its current cluster. Update `pointsPerCluster` for the current cluster label.\n")
  cat("   b. Calculate probabilities for assigning `y[i,]` to each existing cluster `k`:\n")
  cat("      `prob_existing[k] = pointsPerCluster[k] * Likelihood(mdObj, y[i,], clusterParams_k)`\n")
  cat("      where `mdObj` is the mixing distribution object, and `clusterParams_k` are the parameters for cluster `k`.\n")
  cat("   c. Calculate the probability for assigning `y[i,]` to a new cluster:\n")
  cat("      `prob_new = alpha * predictiveArray[i]`\n")
  cat("      where `alpha` is the concentration parameter and `predictiveArray[i]` is the predictive likelihood of `y[i,]` given the base measure.\n")
  cat("   d. Normalize these probabilities (ensure they sum to 1, handle cases where all are zero).\n")
  cat("   e. Sample a new cluster assignment (`newLabel`) for `y[i,]` based on these probabilities.\n")
  cat("   f. Call `ClusterLabelChange()` to update `clusterLabels`, `pointsPerCluster`, `clusterParams`, and `numberClusters` based on `newLabel` and the `currentLabel`. This function handles creating a new cluster if `newLabel` indicates a new cluster, or removing the old cluster if it becomes empty.\n\n")

  cat("## Algorithmic Steps (Non-conjugate Case - `ClusterComponentUpdate.nonconjugate`)\n\n")
  cat("Refers to Neal's Algorithm 8.\n")
  cat("1. For each data point `i` from 1 to `n`:\n")
  cat("   a. Temporarily remove data point `y[i,]` from its current cluster. Update `pointsPerCluster`.\n")
  cat("   b. Calculate probabilities for assigning `y[i,]` to each existing cluster `k`:\n")
  cat("      `prob_existing[k] = pointsPerCluster[k] * Likelihood(mdObj, y[i,], clusterParams_k)`\n")
  cat("   c. If the current cluster becomes empty after removing `y[i,]`, draw `m-1` auxiliary parameters from the prior `PriorDraw(mdObj, m-1)`. Combine these with the parameters of the (now empty) current cluster to form `m` sets of auxiliary parameters `aux`.\n")
  cat("   d. If the current cluster is not empty, draw `m` auxiliary parameters `aux` from `PriorDraw(mdObj, m)`.\n")
  cat("   e. Calculate probabilities for assigning `y[i,]` to new clusters using these auxiliary parameters:\n")
  cat("      `prob_new_aux[j] = (alpha/m) * Likelihood(mdObj, y[i,], aux_j)` for `j` from 1 to `m`.\n")
  cat("   f. Combine `prob_existing` and `prob_new_aux`. Normalize these probabilities.\n")
  cat("   g. Sample a new cluster assignment (`newLabel`) for `y[i,]`.\n")
  cat("   h. Call `ClusterLabelChange()` to update state. If `newLabel` corresponds to an auxiliary parameter, that auxiliary parameter becomes the parameter for the newly formed/assigned cluster.\n\n")

  cat("## Algorithmic Steps (Hierarchical Case - `ClusterComponentUpdate.hierarchical`)\n\n")
  cat("1. Iterates through each individual Dirichlet Process object (`indDP`) in the hierarchical structure.\n")
  cat("2. For each `indDP`, calls the appropriate `ClusterComponentUpdate` method (conjugate or non-conjugate).\n")
  cat("3. After updating components for an `indDP`, calls `DuplicateClusterRemove()` to merge any clusters that may have become identical during the update process.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `dpObj$data`: The input data (matrix).\n")
  cat("- `dpObj$n`: Number of data points.\n")
  cat("- `dpObj$alpha`: Concentration parameter of the DP.\n")
  cat("- `dpObj$clusterLabels`: Vector storing cluster assignment for each data point.\n")
  cat("- `dpObj$clusterParameters`: List of parameters for each cluster.\n")
  cat("- `dpObj$numberClusters`: Current number of active clusters.\n")
  cat("- `dpObj$pointsPerCluster`: Vector storing the count of data points in each cluster.\n")
  cat("- `dpObj$mixingDistribution`: Object containing functions like `Likelihood()`, `PriorDraw()`, `Predictive()`.\n")
  cat("- `dpObj$predictiveArray` (conjugate only): Stores precomputed predictive likelihoods.\n")
  cat("- `dpObj$m` (non-conjugate only): Number of auxiliary parameters.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Sequential Loop:** The primary loop iterates through each data point sequentially, making parallelization for this part challenging.\n")
  cat("2.  **Likelihood Computations:** `Likelihood()` is called multiple times within the loop (for each existing cluster and for new/auxiliary clusters). Efficient likelihood calculation is crucial.\n")
  cat("3.  **`ClusterLabelChange()` Overhead:** This function can involve re-allocating and resizing `clusterParameters` and `pointsPerCluster` if clusters are created or deleted, which can be costly.\n")
  cat("4.  **`PriorDraw()` (non-conjugate):** Drawing `m` auxiliary parameters in each iteration adds computational load.\n")
  cat("5.  **S3 Dispatch:** Minor overhead from S3 method dispatch, but likely less significant than the above points.\n\n")

  cat("## Dependencies\n\n")
  cat("- `Likelihood()` method for the specific mixing distribution.\n")
  cat("- `Predictive()` method (conjugate case) for the mixing distribution.\n")
  cat("- `PriorDraw()` method (non-conjugate case) for the mixing distribution.\n")
  cat("- `ClusterLabelChange()` function (handles the logic of updating cluster assignments, parameters, and counts after a new label is chosen).\n")
  cat("- `DuplicateClusterRemove()` function (hierarchical case).\n\n")
  sink()
}

analyze_cluster_parameter_update <- function() {
  sink("analysis/code_documentation/cluster_parameter_update.md")

  cat("# ClusterParameterUpdate Algorithm Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("This function updates the parameters for each existing cluster based on the data points currently assigned to it. It's a crucial step in the MCMC fitting process of a Dirichlet Process Mixture Model (DPMM).\n\n")

  cat("## Implementation Structure\n\n")
  cat("The function is implemented with S3 dispatch based on the `dpObj` class (specifically, its conjugacy status):\n")
  cat("- `ClusterParameterUpdate()`: Generic function.\n")
  cat("- `ClusterParameterUpdate.conjugate()`: Handles DPMMs with conjugate base measures.\n")
  cat("- `ClusterParameterUpdate.nonconjugate()`: Handles DPMMs with non-conjugate base measures.\n\n")
  cat("Note: There isn't a specific `ClusterParameterUpdate.hierarchical()` method explicitly shown in `cluster_parameter_update.R`. Hierarchical parameter updates are likely handled within `GlobalParameterUpdate` and updates to individual DPs within the hierarchy would use their respective conjugate/non-conjugate methods. The proposal mentions `GlobalParameterUpdate` for hierarchical models (Section 8.1 of priyanshuTiwari_dirichletprocess.pdf).\n\n")

  cat("## Algorithmic Steps (Conjugate Case - `ClusterParameterUpdate.conjugate`)\n\n")
  cat("1. For each existing cluster `i` from 1 to `numLabels` (current number of clusters):\n")
  cat("   a. Identify all data points `pts` assigned to cluster `i` (i.e., where `clusterLabels == i`).\n")
  cat("   b. Draw new parameters for cluster `i` from the posterior distribution: `post_draw = PosteriorDraw(mdobj, pts)`.\n")
  cat("      `mdobj` is the mixing distribution object. `PosteriorDraw` leverages conjugacy to sample directly from the analytical posterior given the data `pts` assigned to the cluster.\n")
  cat("   c. Update `clusterParams` for cluster `i` with these `post_draw` values.\n\n")

  cat("## Algorithmic Steps (Non-conjugate Case - `ClusterParameterUpdate.nonconjugate`)\n\n")
  cat("1. For each existing cluster `i` from 1 to `numLabels`:\n")
  cat("   a. Identify all data points `pts` assigned to cluster `i`.\n")
  cat("   b. Get the current parameters for cluster `i` to use as a starting position for the Metropolis-Hastings sampler: `start_pos`.\n")
  cat("   c. Draw `mhDraws` samples from the posterior distribution using Metropolis-Hastings: `parameter_samples = PosteriorDraw(mdobj, pts, mhDraws, start_pos = start_pos)`.\n")
  cat("      `PosteriorDraw` for non-conjugate cases internally uses `MetropolisHastings()` which involves proposing new parameters (`MhParameterProposal()`), calculating likelihoods (`Likelihood()`), prior densities (`PriorDensity()`), and acceptance probabilities.\n")
  cat("   d. Update `clusterParams` for cluster `i` with the last sample from `parameter_samples`.\n")
  cat("   e. (Optionally) Calculate and store the acceptance ratio for the MH step for diagnostics.\n\n")

  cat("The internal function `cluster_parameter_update` (not exported, and seemingly a conceptual or alternative implementation) also outlines a similar process where it iterates unique clusters and calls `PosteriorDraw` for data belonging to each.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `dpObj$data`: The input data.\n")
  cat("- `dpObj$clusterLabels`: Vector of cluster assignments.\n")
  cat("- `dpObj$clusterParameters`: List of parameters for each cluster (this is what gets updated).\n")
  cat("- `dpObj$numberClusters`: Current number of active clusters.\n")
  cat("- `dpObj$mixingDistribution`: Object containing `PosteriorDraw()`, `PriorDraw()`, `Likelihood()`, `PriorDensity()`, `MhParameterProposal()` methods.\n")
  cat("- `dpObj$mhDraws` (non-conjugate only): Number of Metropolis-Hastings iterations for parameter sampling.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Looping over Clusters:** The algorithm iterates through each active cluster.\n")
  cat("2.  **Data Subsetting:** For each cluster, it subsets the data points belonging to it.\n")
  cat("3.  **`PosteriorDraw()` Cost:**\n")
  cat("    -   **Conjugate:** Relatively efficient as it involves direct sampling from a known posterior distribution. The complexity depends on calculating posterior parameters from data in the cluster.\n")
  cat("    -   **Non-conjugate:** Can be significantly more expensive due to the iterative nature of `MetropolisHastings()`. This involves `mhDraws` iterations of proposing parameters, and repeatedly calculating likelihoods and prior densities.\n")
  cat("4.  **Likelihood and Prior Density Calls (Non-conjugate):** Inside `MetropolisHastings`, `Likelihood()` and `PriorDensity()` are called for each of the `mhDraws` iterations per cluster.\n\n")

  cat("## Dependencies\n\n")
  cat("- `PosteriorDraw()` method for the specific mixing distribution (which itself depends on other methods for non-conjugate cases, such as `MetropolisHastings`, `MhParameterProposal`, `Likelihood`, `PriorDensity`).\n\n")
  sink()
}

analyze_update_alpha <- function() {
  sink("analysis/code_documentation/update_alpha.md")

  cat("# UpdateAlpha Algorithm Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("This function updates the concentration parameter `alpha` of the Dirichlet Process. A well-chosen `alpha` is crucial for determining the number of clusters in the DPMM. The update is typically based on the current number of clusters and the total number of data points, following methods like those described by West (1992).\n\n")

  cat("## Implementation Structure\n\n")
  cat("The function is implemented with S3 dispatch:\n")
  cat("- `UpdateAlpha()`: Generic function.\n")
  cat("- `UpdateAlpha.default()`: Handles the standard DP case.\n")
  cat("- `UpdateAlpha.hierarchical()`: Handles hierarchical DPs by iterating through individual DPs and updating their respective alpha values, and also updating their `mixingDistribution$alpha` to match.\n\n")

  cat("## Algorithmic Steps (`UpdateAlpha.default` via `update_concentration`)\n\n")
  cat("The core logic is in the internal function `update_concentration(oldParam, n, nParams, priorParameters)` where `oldParam` is the current alpha, `n` is the total number of data points, `nParams` is `dpObj$numberClusters`, and `priorParameters` are `dpObj$alphaPriorParameters` (shape `a` and rate `b` for a Gamma prior on alpha).\n")
  cat("The method is based on Escobar and West (1995), introducing an auxiliary variable.\n")
  cat("1. Draw an auxiliary variable `x` from a Beta distribution: `x ~ Beta(alpha + 1, n)`.\n")
  cat("2. Define weights for a mixture of two Gamma distributions:\n")
  cat("   `pi1_numerator = priorParameters[1] + nParams - 1` (shape_prior + K - 1)\n")
  cat("   `pi1_denominator_term = n * (priorParameters[2] - log(x))` (n * (rate_prior - log(x)))\n")
  cat("   `pi1 = pi1_numerator / (pi1_numerator + pi1_denominator_term)`\n")
  cat("3. Define posterior parameters for the Gamma distributions:\n")
  cat("   `postParams1_shape_candidate1 = priorParameters[1] + nParams`\n")
  cat("   `postParams1_shape_candidate2 = priorParameters[1] + nParams - 1`\n")
  cat("   `postParams2_rate = priorParameters[2] - log(x)`\n")
  cat("4. With probability `pi1`, the new alpha is drawn from `Gamma(postParams1_shape_candidate1, postParams2_rate)`.\n")
  cat("5. Otherwise (with probability `1-pi1`), the new alpha is drawn from `Gamma(postParams1_shape_candidate2, postParams2_rate)`.\n")
  cat("   (The R code slightly simplifies this by adjusting `postParams1_shape_candidate1` if `runif(1) > pi1`).\n")
  cat("6. The drawn value becomes the new `dpObj$alpha`.\n\n")

  cat("## Algorithmic Steps (`UpdateAlpha.hierarchical`)\n\n")
  cat("1. For each individual Dirichlet Process (`indDP`) within the `dpObjList`:\n")
  cat("   a. Call `UpdateAlpha(dpObj$indDP[[i]])` (which will dispatch to `UpdateAlpha.default` for that individual DP).\n")
  cat("   b. Update `dpObj$indDP[[i]]$mixingDistribution$alpha = dpObj$indDP[[i]]$alpha` to ensure consistency.\n\n")
  cat("A similar function `UpdateGamma` exists for the global concentration parameter in hierarchical models, using the same `update_concentration` logic but with global counts of tables and parameters.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `dpObj$alpha`: The current concentration parameter (updated).\n")
  cat("- `dpObj$n`: Total number of data points.\n")
  cat("- `dpObj$numberClusters`: Current number of active clusters.\n")
  cat("- `dpObj$alphaPriorParameters`: A vector `c(shape, rate)` for the Gamma prior on alpha.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Random Number Generation:** Involves drawing from Beta and Gamma distributions, which are generally efficient.\n")
  cat("2.  **Logarithm and Basic Arithmetic:** Calculations are primarily arithmetic and include one log computation.\n")
  cat("3.  **Overall Cost:** This step is typically much faster than `ClusterComponentUpdate` or `ClusterParameterUpdate` as it does not loop over data points or clusters, only performing a few calculations and random draws.\n\n")

  cat("## Dependencies\n\n")
  cat("- `rbeta()`: For drawing from the Beta distribution.\n")
  cat("- `rgamma()`: For drawing from the Gamma distribution.\n")
  cat("- `runif()`: For deciding which component of the mixture Gamma to draw from.\n\n")

  sink()
}

analyze_likelihood_functions <- function() {
  sink("analysis/code_documentation/likelihood_functions.md")

  cat("# Likelihood Functions (`Likelihood.*`) Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("The `Likelihood()` functions are S3 generic methods responsible for calculating the likelihood of one or more data points `x` given a specific set of parameters `theta` for a particular mixing distribution (kernel) defined by `mdObj`. These are fundamental for the DPMM fitting process, used in assigning points to clusters and in updating cluster parameters (especially in non-conjugate cases).\n\n")

  cat("## Implementation Structure\n\n")
  cat("S3 dispatch is used based on the class of the mixing distribution object (`mdObj`). Specific methods exist for each supported kernel:\n")
  cat("- `Likelihood.normal()`: For Gaussian mixture kernel (unknown mean and variance).\n")
  cat("- `Likelihood.normalFixedVariance()`: For Gaussian mixture kernel with fixed variance.\n")
  cat("- `Likelihood.beta()`: For Beta mixture kernel (parameterized by mean and precision, on `[0, maxT]`).\n")
  cat("- `Likelihood.beta2()`: For Beta mixture kernel (alternative parameterization/prior, also on `[0, maxT]`).\n")
  cat("- `Likelihood.weibull()`: For Weibull mixture kernel.\n")
  cat("- `Likelihood.exponential()`: For Exponential mixture kernel.\n")
  cat("- `Likelihood.mvnormal()`: For Multivariate Normal mixture kernel (conjugate Normal-Wishart prior structure).\n")
  cat("- `Likelihood.mvnormal2()`: For Multivariate Normal mixture kernel (semi-conjugate, independent priors on mean and covariance).\n")
  cat("- Custom implementations for user-defined kernels would follow this pattern.\n\n")

  cat("## Algorithmic Steps (General)\n\n")
  cat("For a given data point(s) `x` and parameter set `theta`:\n")
  cat("1. Extract the relevant parameters from the `theta` list (e.g., mean and standard deviation for Normal, shape and scale for Weibull, etc.). The `theta` argument is often a list where each element is an array, allowing for vectorized calculations if multiple parameter sets are provided (e.g., for different clusters or auxiliary parameters).\n")
  cat("2. Handle the dimensionality: `theta` parameters are often structured as 3D arrays `(param_dim1, param_dim2, num_param_sets)`. The function needs to correctly extract parameters for each set if `x` is being evaluated against multiple `theta`.\n")
  cat("3. For each data point in `x` (if `x` is a vector or matrix of multiple points) and for each parameter set in `theta`:\n")
  cat("   a. Call the corresponding base R density function (e.g., `dnorm()`, `dbeta()`, `dweibull()`, `mvtnorm::dmvnorm()`).\n")
  cat("   b. Adjustments for parameterization might be needed. For example, `Likelihood.beta` converts mean/precision (mu, tau) to shape1/shape2 (a,b) and scales by `1/maxT` because `dbeta` in R is for the standard [0,1] Beta.\n")
  cat("4. Return a numeric vector or matrix of likelihood values. If multiple parameter sets in `theta` are provided, the output is typically a vector where each element corresponds to the likelihood of `x` under one parameter set (e.g., in `Likelihood.mvnormal`).\n\n")

  cat("## Specific Kernel Details\n\n")
  cat("### `Likelihood.normal` (normal_inverse_gamma.R)\n")
  cat("- Parameters in `theta`: `[[1]]` is mean (`mu`), `[[2]]` is standard deviation (`sigma`).\n")
  cat("- Uses `dnorm(x, theta[[1]], theta[[2]])`.\n\n")

  cat("### `Likelihood.beta` (beta_uniform_gamma.R)\n")
  cat("- Parameters in `theta`: `[[1]]` is mean (`mu`), `[[2]]` is precision/scale (`tau` or `nu`).\n")
  cat("- `maxT` is retrieved from `mdObj$maxT`.\n")
  cat("- Converts (mu, tau) to standard Beta parameters (a,b): `a = (mu * tau)/maxT`, `b = (1 - mu/maxT) * tau`.\n")
  cat("- Calculates likelihood as `1/maxT * dbeta(x/maxT, a, b)` to account for the `[0, maxT]` support.\n\n")

  cat("### `Likelihood.weibull` (weibull_uniform_gamma.R)\n")
  cat("- Parameters in `theta`: `[[1]]` is shape (`alpha`), `[[2]]` is scale (`lambda`). The code uses $ \\lambda^{-1} \\alpha x^{\\alpha-1} \\exp(-\\lambda^{-1} x^\\alpha) $ which means $\\lambda$ is the scale parameter `b`.\n") # Escaped LaTeX
  cat("- Handles `Inf` values for `lambda` (scale) by setting likelihood to 0.\n")
  cat("- Sets likelihood to 0 for `x < 0`.\n\n")

  cat("### `Likelihood.mvnormal` & `Likelihood.mvnormal2` (mvnormal_normal_wishart.R, mvnormal_semi_conjugate.R)\n")
  cat("- Parameters in `theta`: `[[1]]` (or `theta$mu`) is the mean vector, `[[2]]` (or `theta$sig`) is the covariance matrix.\n")
  cat("- Uses `mvtnorm::dmvnorm(x, mean, sigma)`.\n")
  cat("- Loops or `vapply`s through multiple parameter sets if provided in `theta`.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `mdObj`: The mixing distribution object, which might contain fixed parameters (e.g., `mdObj$sigma` in `Likelihood.normalFixedVariance`, `mdObj$maxT` in `Likelihood.beta`).\n")
  cat("- `x`: The data point(s) (numeric vector or matrix).\n")
  cat("- `theta`: A list of parameter arrays. Each element of the list corresponds to a parameter of the distribution (e.g., `theta[[1]]` for $\\mu$, `theta[[2]]` for $\\sigma^2$). The third dimension of these arrays typically allows for evaluating likelihood against multiple parameter sets.\n\n") # Escaped LaTeX

  cat("## Performance Considerations\n\n")
  cat("1.  **Vectorization**: Efficient implementations vectorize calculations over `x` if `theta` contains a single parameter set, or over `theta` if `x` is a single data point evaluated against multiple parameter sets. The use of `vapply` in some mvnormal cases suggests iteration over parameter sets.\n")
  cat("2.  **Base R Efficiency**: Relies on the efficiency of base R's density functions (`dnorm`, `dbeta`, etc.) or `mvtnorm::dmvnorm()`.\n")
  cat("3.  **Parameter Transformation**: For some kernels (e.g., Beta), parameters might need transformation before calling the base R density function, adding slight overhead.\n")
  cat("4.  **Repeated Calls**: These functions are called very frequently. Even small inefficiencies can compound.\n\n")

  cat("## Dependencies\n\n")
  cat("- Base R density functions (e.g., `dnorm`, `dbeta`, `dweibull`, `dexp`).\n")
  cat("- `mvtnorm::dmvnorm` for multivariate normal distributions.\n\n")

  sink()
}

analyze_prior_draw_functions <- function() {
  sink("analysis/code_documentation/prior_draw_functions.md")

  cat("# PriorDraw Functions (`PriorDraw.*`) Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("The `PriorDraw()` S3 generic methods are responsible for drawing `n` samples of parameters `theta` from the base measure `G_0` associated with a specific mixing distribution `mdObj`. These draws are used to parameterize new clusters when they are formed, and to generate auxiliary parameters in Neal's Algorithm 8 for non-conjugate mixtures.\n\n")

  cat("## Implementation Structure\n\n")
  cat("S3 dispatch is used based on the class of `mdObj`. Each mixing distribution has its own method defining how to sample from its specific base measure `G_0`:\n")
  cat("- `PriorDraw.normal()`: Normal-Inverse Gamma base measure for Gaussian kernel.\n")
  cat("- `PriorDraw.normalFixedVariance()`: Normal base measure for mean (variance is fixed).\n")
  cat("- `PriorDraw.beta()`: Uniform for mean `mu`, Inverse-Gamma for scale `nu`.\n")
  cat("- `PriorDraw.beta2()`: Uniform for mean `mu`, Pareto for scale `nu`.\n")
  cat("- `PriorDraw.weibull()`: Uniform for shape `a`, Inverse-Gamma for scale `lambda` (used as `b`).\n")
  cat("- `PriorDraw.exponential()`: Gamma base measure for the rate parameter.\n")
  cat("- `PriorDraw.mvnormal()`: Normal for mean `mu`, Wishart for precision matrix `sig` (likely covariance in practice).\n") # Adjusted based on typical use of rWishart
  cat("- `PriorDraw.mvnormal2()`: Normal for mean `mu`, Inverse-Wishart for covariance matrix `sig`.\n")
  cat("- `PriorDraw.hierarchical()`: Samples from the global discrete measure $G_0$ which itself is a draw from a DP. It samples an index based on global stick-breaking weights `pi_k` and returns the corresponding global parameter `theta_k`.\n\n")

  cat("## Algorithmic Steps (General)\n\n")
  cat("1. Retrieve prior parameters from `mdObj$priorParameters`.\n")
  cat("2. For each of the `n` samples requested:\n")
  cat("   a. Draw each component of the parameter vector `theta` from its specified prior distribution using base R's random number generators (e.g., `rnorm()`, `rgamma()`, `runif()`, `rWishart()`).\n")
  cat("   b. The parameters drawn might be independent or dependent (e.g., in Normal-Inverse Gamma, the variance is drawn first, then the mean conditional on the variance).\n")
  cat("3. Structure the `n` samples into a list of arrays. Each list element corresponds to a parameter (e.g., `theta[[1]]` for $\\mu$, `theta[[2]]` for $\\sigma$), and each array has dimensions appropriate for the parameter (e.g., `1x1xn` for univariate scalar parameters, `dx1xn` for mean vectors, `dxdxn` for covariance matrices).\n\n") # Escaped LaTeX

  cat("## Specific Kernel Details\n\n")
  cat("### `PriorDraw.normal` (normal_inverse_gamma.R)\n")
  cat("- Draws `lambda` (precision) from `Gamma(alpha0, beta0)`.\n")
  cat("- Draws `mu` (mean) from `Normal(mu0, (kappa0 * lambda)^-0.5)`.\n")
  cat("- Returns `list(array(mu), array(sqrt(1/lambda)))`.\n\n")

  cat("### `PriorDraw.beta` (beta_uniform_gamma.R)\n")
  cat("- Draws `mu` from `Uniform(0, mdObj$maxT)`.\n")
  cat("- Draws `nu` (scale) from `1 / Gamma(alpha0, beta0)` (i.e., Inverse-Gamma).\n\n")

  cat("### `PriorDraw.weibull` (weibull_uniform_gamma.R)\n")
  cat("- Draws shape `a` from `Uniform(0, priorParameters[1])` (phi).\n")
  cat("- Draws scale `lambda` (representing `b`) from `1 / Gamma(priorParameters[2], priorParameters[3])` (alpha0, beta0 for Inv-Gamma).\n\n")

  cat("### `PriorDraw.mvnormal` (mvnormal_normal_wishart.R)\n")
  cat("- Draws covariance matrix `sig` from `Wishart(nu, Lambda)` (where `Lambda` is the scale matrix for `rWishart`).\n") # Clarified based on rWishart
  cat("- Draws mean `mu` from `Normal(mu0, solve(sig * kappa0))` (multivariate normal).\n\n")

  cat("### `PriorDraw.hierarchical` (mixing_distribution_prior_draw.R)\n")
  cat("- Uses probabilities `mdObj$pi_k` (global stick-breaking weights for $G_0$).\n")
  cat("- Samples `n` indices `ind` based on these probabilities.\n")
  cat("- Returns `lapply(mdObj$theta_k, function(x) x[,, ind, drop = FALSE])`, effectively selecting `n` parameters from the existing discrete global atoms `theta_k`.\n\n")


  cat("## Key Data Structures Used\n\n")
  cat("- `mdObj$priorParameters`: Contains the hyperparameters for the base measure `G_0`.\n")
  cat("- `n`: The number of parameter sets to draw.\n")
  cat("- Output: A list of arrays, where each array holds `n` samples of a specific parameter.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Efficiency of RNGs**: Performance depends on the efficiency of R's underlying random number generators for distributions like Gamma, Normal, Uniform, Wishart.\n")
  cat("2.  **Looping for `n` samples**: Each of the `n` samples is typically drawn independently.\n")
  cat("3.  **Complexity of Prior**: More complex prior structures (e.g., hierarchical dependencies between parameters) will take longer to sample.\n")
  cat("4.  **`PriorDraw.hierarchical`**: Depends on sampling from a discrete distribution, which is generally efficient.\n\n")

  cat("## Dependencies\n\n")
  cat("- Base R random number generation functions (`rgamma`, `rnorm`, `runif`, etc.).\n")
  cat("- `rWishart` from the `stats` package for multivariate normal base measures.\n\n")

  sink()
}

analyze_posterior_draw_functions <- function() {
  sink("analysis/code_documentation/posterior_draw_functions.md")

  cat("# PosteriorDraw Functions (`PosteriorDraw.*`) Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("The `PosteriorDraw()` S3 generic methods are responsible for drawing `n` samples of parameters `theta` from the posterior distribution $p(\\theta | x, G_0)$, given the data `x` assigned to a specific cluster and the base measure `G_0` (via `mdObj`). This is central to the `ClusterParameterUpdate` step.\n\n") # Escaped LaTeX

  cat("## Implementation Structure\n\n")
  cat("S3 dispatch is used based on the class of `mdObj`, distinguishing mainly between conjugate and non-conjugate mixing distributions:\n")
  cat("- `PosteriorDraw.conjugate` (implicit, specific methods for each conjugate kernel):\n")
  cat("  - `PosteriorDraw.normal()` (normal_inverse_gamma.R): Samples from Normal-Inverse Gamma posterior.\n")
  cat("  - `PosteriorDraw.normalFixedVariance()` (normal_fixed_variance.R): Samples from Normal posterior for the mean.\n")
  cat("  - `PosteriorDraw.exponential()` (exponential_gamma.R): Samples from Gamma posterior for the rate.\n")
  cat("  - `PosteriorDraw.mvnormal()` (mvnormal_normal_wishart.R): Samples from Normal-Wishart posterior.\n")
  cat("- `PosteriorDraw.nonconjugate` (mixing_distribution_posterior_draw.R): Generic for non-conjugate. This method typically wraps `MetropolisHastings()`.\n")
  cat("  - Specific non-conjugate kernels like `PosteriorDraw.beta()` (beta_uniform_gamma.R), `PosteriorDraw.weibull()` (weibull_uniform_gamma.R), `PosteriorDraw.mvnormal2()` (mvnormal_semi_conjugate.R) implement this by calling `MetropolisHastings` or having their own MH-like loop.\n\n")

  cat("## Algorithmic Steps (Conjugate Cases)\n\n")
  cat("1. Calculate posterior hyperparameters: The function first calls `PosteriorParameters(mdObj, x)` to compute the parameters of the analytical posterior distribution based on the prior parameters from `mdObj` and the data `x` in the current cluster.\n")
  cat("2. Draw `n` samples from this analytical posterior distribution using base R random number generators.\n")
  cat("   - E.g., for `PosteriorDraw.normal()`: Draw `lambda` (precision) from `Gamma(alpha_n, beta_n)`, then `mu` from `Normal(mu_n, (kappa_n * lambda)^-0.5)`.\n")
  cat("3. Structure the samples into the standard list of arrays format.\n\n")

  cat("## Algorithmic Steps (Non-conjugate Cases - General via `PosteriorDraw.nonconjugate`)\n\n")
  cat("1. Determine starting position for Metropolis-Hastings: \n")
  cat("   - If `start_pos` is provided via `...`, use it.\n")
  cat("   - Else, call `PenalisedLikelihood(mdObj, x)` to find a good starting point (e.g., mode of a penalized likelihood) or fall back to `PriorDraw(mdObj, 1)` if `PenalisedLikelihood` is not implemented for the specific `mdObj`.\n")
  cat("2. Call `MetropolisHastings(mdObj, x, start_pos, no_draws = n)` to obtain `n` samples from the posterior via MCMC.\n")
  cat("3. Format the output samples from `MetropolisHastings` into the standard list of arrays.\n\n")

  cat("### Specific Non-conjugate Implementations (e.g., `PosteriorDraw.weibull`):\n")
  cat("- May have custom logic within the MH loop, for example, if one parameter can be updated via Gibbs sampling conditional on the other MH-sampled parameter (as in Weibull where `lambda` is sampled from its conditional posterior given `alpha`).\n") # Escaped LaTeX
  cat("- `PosteriorDraw.weibull` explicitly samples `lambda` from its conditional posterior (an Inverse Gamma derived from the Gamma on $1/\\lambda$) given the current `alpha` proposed by `MhParameterProposal.weibull`.\n\n") # Escaped LaTeX

  cat("## Key Data Structures Used\n\n")
  cat("- `mdObj`: Mixing distribution object, containing prior parameters and methods like `Likelihood`, `PriorDensity`, `MhParameterProposal`.\n")
  cat("- `x`: Data points assigned to the current cluster.\n")
  cat("- `n`: Number of posterior samples to draw.\n")
  cat("- `...` (often `start_pos`): Starting parameters for non-conjugate MCMC.\n")
  cat("- Output: A list of arrays, similar to `PriorDraw` output.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Conjugate Models**: Performance is dictated by the cost of calculating posterior hyperparameters (usually involves summing over data `x`) and then drawing `n` samples. Generally efficient.\n")
  cat("2.  **Non-conjugate Models**: Significantly more computationally intensive.\n")
  cat("    -   Dominated by the `MetropolisHastings` call, which runs for `n` (often `mhDraws`) iterations.\n")
  cat("    -   Each MH iteration involves: parameter proposal, `Likelihood()` call (summed over data `x`), `PriorDensity()` call, and acceptance step.\n")
  cat("    -   The cost of `Likelihood(mdObj, x, theta)` is crucial here as `x` can be many data points.\n")
  cat("    -   `PenalisedLikelihood()` (if used for `start_pos`) can involve an optimization routine which adds to the setup cost.\n\n")

  cat("## Dependencies\n\n")
  cat("- For conjugate cases: `PosteriorParameters()` method, base R RNGs.\n")
  cat("- For non-conjugate cases: `MetropolisHastings()` generic (and its specific methods), `MhParameterProposal()`, `Likelihood()`, `PriorDensity()`, and potentially `PenalisedLikelihood()` or `PriorDraw()` for starting positions.\n\n")

  sink()
}

analyze_metropolis_hastings_components <- function() {
  sink("analysis/code_documentation/metropolis_hastings_components.md")

  cat("# Metropolis-Hastings Components Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("These components form the core of posterior sampling for non-conjugate Dirichlet Process Mixture Models. They are used within `PosteriorDraw()` for non-conjugate kernels to sample from the posterior distribution of cluster parameters $p(\\theta | x, G_0)$, where `x` is the data assigned to a particular cluster.\n\n") # Escaped LaTeX
  cat("The main functions involved are:\n")
  cat("- `MetropolisHastings()` (S3 generic with methods like `MetropolisHastings.default` and `MetropolisHastings.weibull`)\n")
  cat("- `MhParameterProposal()` (S3 generic, e.g., `MhParameterProposal.beta`, `MhParameterProposal.weibull`)\n")
  cat("- `PriorDensity()` (S3 generic, e.g., `PriorDensity.beta`, `PriorDensity.weibull`)\n")
  cat("- `Likelihood()` (S3 generic, already analyzed separately, but crucial here)\n\n")

  cat("## `MetropolisHastings.*` Methods\n\n")
  cat("### Algorithm Overview (`MetropolisHastings.default` from `metropolis_hastings.R`)\n")
  cat("1. Initialize: Set `parameter_samples` array. The first sample is the `start_pos`.\n")
  cat("2. Calculate initial log-posterior components: `old_prior = log(PriorDensity(mixingDistribution, old_param))` and `old_Likelihood = sum(log(Likelihood(mixingDistribution, x, old_param)))`.\n")
  cat("3. For `i` from 1 to `no_draws - 1`:\n")
  cat("   a. Propose new parameters: `prop_param = MhParameterProposal(mixingDistribution, old_param)`.\n")
  cat("   b. Calculate log-posterior components for the proposal: `new_prior = log(PriorDensity(mixingDistribution, prop_param))` and `new_Likelihood = sum(log(Likelihood(mixingDistribution, x, prop_param)))`.\n")
  cat("   c. Calculate acceptance probability `accept_prob = min(1, exp(new_prior + new_Likelihood - old_prior - old_Likelihood))`.\n")
  cat("      - Handle `NA` or invalid `accept_prob` by setting it to 0.\n")
  cat("   d. If `runif(1) < accept_prob` (accept proposal):\n")
  cat("      - `sampled_param = prop_param`\n")
  cat("      - `old_Likelihood = new_Likelihood`\n")
  cat("      - `old_prior = new_prior`\n")
  cat("   e. Else (reject proposal):\n")
  cat("      - `sampled_param = old_param`\n")
  cat("   f. Store `sampled_param` in `parameter_samples`.\n")
  cat("   g. Update `old_param = sampled_param`.\n")
  cat("4. Return `list(parameter_samples, accept_ratio)`.\n\n")

  cat("### `MetropolisHastings.weibull` (metropolis_hastings.R)\n")
  cat("- This is a specialized version that incorporates a Gibbs step for one of the parameters.\n")
  cat("- After proposing the shape parameter `alpha` (from `MhParameterProposal.weibull`), it directly samples the scale parameter `lambda` from its full conditional posterior `Gamma(length(x) + priorParameters[2], sum(x^alpha_prop) + priorParameters[3])`, (actually sampling $1/\\lambda$ then inverting to get $\\lambda$).\n") # Escaped LaTeX
  cat("- The acceptance probability is then based on the joint proposal (MH for `alpha`, Gibbs for `lambda`).\n\n") # Escaped LaTeX

  cat("## `MhParameterProposal.*` Methods\n\n")
  cat("### Algorithm Overview\n")
  cat("These S3 methods define how a new candidate parameter set `new_params` is proposed based on the `old_params` and `mdObj$mhStepSize`.\n")
  cat("- Typically involve adding random noise (e.g., from a Normal distribution scaled by `mhStepSize`) to the `old_params`.\n")
  cat("- Must respect parameter constraints (e.g., positivity for scale/variance parameters, bounds for probabilities).\n\n")
  cat("### Specific Implementations\n")
  cat("- `MhParameterProposal.beta` (beta_uniform_gamma.R): Adds Normal noise to `mu` (mean), clips to `[0, maxT]`. Adds Normal noise to `nu` (scale/precision) and takes absolute value to ensure positivity.\n")
  cat("- `MhParameterProposal.weibull` (weibull_uniform_gamma.R): Adds Normal noise to shape `alpha` and takes absolute value.\n") # Escaped LaTeX
  cat("- `MhParameterProposal.beta2` (beta_uniform_pareto.R): Similar to `.beta` but for potentially different step sizes or proposal logic if needed.\n\n")

  cat("## `PriorDensity.*` Methods\n\n")
  cat("### Algorithm Overview\n")
  cat("These S3 methods calculate the log-prior density $\\log p(\\theta)$ for a given parameter set `theta` based on the `mdObj$priorParameters`.\n\n") # Escaped LaTeX
  cat("### Specific Implementations\n")
  cat("- `PriorDensity.beta` (beta_uniform_gamma.R): `dunif(mu, 0, maxT) * dgamma(1/nu, shape, rate)` (product of densities for mu and nu).\n")
  cat("- `PriorDensity.weibull` (weibull_uniform_gamma.R): `dunif(shape_a, 0, phi) * dgamma(1/scale_b, alpha0, beta0)`. (The R code only shows `dunif(theta[[1]], 0, priorParameters[1])` suggesting the scale part might be integrated out or handled differently in the Weibull MH or it's a simplified prior density for MH purposes, which is common if one parameter has a direct sampler).\n")
  cat("- `PriorDensity.beta2` (beta_uniform_pareto.R): `dunif(mu, 0, maxT) * dpareto(nu, muLim, gamma_shape)`.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `mixingDistribution` (`mdObj`): Contains prior parameters, `mhStepSize`, and kernel type.\n")
  cat("- `x`: Data points for likelihood calculation.\n")
  cat("- `start_pos`: Initial parameter values for the MH chain.\n")
  cat("- `no_draws`: Number of MH iterations.\n")
  cat("- `old_param`, `prop_param`, `sampled_param`: Parameter lists/vectors.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Number of Draws (`no_draws`):** The main loop runs this many times. This is often set by `dpObj$mhDraws`.\n")
  cat("2.  **`Likelihood()` Call:** This is the most expensive part within each MH iteration, especially if the number of data points `x` in the cluster is large.\n")
  cat("3.  **`PriorDensity()` Call:** Usually less expensive than `Likelihood()`.\n")
  cat("4.  **`MhParameterProposal()` Call:** Typically very fast (random draws and arithmetic).\n")
  cat("5.  **Log/Exp and Arithmetic:** Standard arithmetic operations.\n")
  cat("6.  **Acceptance Rate:** Low acceptance rates mean many proposals (and their likelihood/prior calculations) are wasted, leading to inefficiency and poor mixing.\n\n")

  cat("## Dependencies\n\n")
  cat("- `Likelihood.*()` methods for the specific kernel.\n")
  cat("- `PriorDensity.*()` methods for the specific kernel.\n")
  cat("- `MhParameterProposal.*()` methods for the specific kernel.\n")
  cat("- Base R random number generators (e.g., `rnorm`, `runif`).\n\n")
  sink()
}

analyze_predictive_likelihood_functions <- function() {
  sink("analysis/code_documentation/predictive_likelihood_functions.md")

  cat("# Predictive Likelihood Functions (`Predictive.*`) Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("The `Predictive()` S3 generic methods calculate the marginal likelihood of new data point(s) `x` given the base measure `G_0` (defined by `mdObj`). This is $p(x | G_0) = \\int k(x | \\theta) dG_0(\\theta)$. These functions are primarily used in the `ClusterComponentUpdate.conjugate()` method to determine the probability of assigning a data point to a new cluster.\n\n") # Escaped LaTeX

  cat("## Implementation Structure\n\n")
  cat("S3 dispatch is used based on the class of the mixing distribution object (`mdObj`). These are typically implemented only for conjugate mixture models, as the integral is often analytically tractable in such cases:\n")
  cat("- `Predictive.normal()` (normal_inverse_gamma.R): For the Normal-Inverse Gamma base measure.\n")
  cat("- `Predictive.normalFixedVariance()` (normal_fixed_variance.R): For the Normal base measure (fixed variance).\n")
  cat("- `Predictive.exponential()` (exponential_gamma.R): For the Gamma base measure.\n")
  cat("- `Predictive.mvnormal()` (mvnormal_normal_wishart.R): For the Normal-Wishart base measure.\n")
  cat("- Non-conjugate kernels generally do not have a simple `Predictive()` method because the integral $\\int k(x | \\theta) dG_0(\\theta)$ is not analytically tractable. In such cases, Neal's Algorithm 8 (using auxiliary parameters) is used instead of relying on this predictive likelihood for new cluster formation.\n\n") # Escaped LaTeX

  cat("## Algorithmic Steps (General for Conjugate Kernels)\n\n")
  cat("1. Retrieve prior parameters from `mdObj$priorParameters`.\n")
  cat("2. For each data point `x[i]` (if `x` is a vector/matrix of multiple points):\n")
  cat("   a. Analytically derive or use known results for the marginal likelihood (predictive distribution) of `x[i]` under the specified `mdObj` (kernel `k` and base measure `G_0`).\n")
  cat("   b. This often involves updating the prior parameters to 'posterior' parameters as if `x[i]` was observed, and then using ratios of normalizing constants or properties of standard distributions.\n")
  cat("      - For example, in `Predictive.normal()`, it calculates the posterior parameters as if `x[i]` was the only data point, then uses a formula involving Gamma functions and ratios of prior/posterior parameters (related to the Student-t predictive distribution).\n")
  cat("      - For `Predictive.exponential()`, it computes $p(x_i) = \\frac{\\Gamma(\\alpha_0+1)}{\\Gamma(\\alpha_0)} \\frac{\\beta_0^{\\alpha_0}}{(\\beta_0+x_i)^{\\alpha_0+1}}$ (if $x_i$ is a single observation and length(x[i])=1, and sum(x[i])=x[i]). The code `alphaPost <- priorParameters[1] + length(x[i]); betaPost <- priorParameters[2] + sum(x[i])` suggests it's calculating for a set of x[i] as if they formed a single new observation, which is correct for the predictive likelihood of a new data point being integrated over the prior.\n") # Escaped LaTeX
  cat("3. Store the calculated predictive likelihood for each `x[i]` in `predictiveArray`.\n")
  cat("4. Return `predictiveArray`.\n\n")

  cat("## Specific Kernel Details\n\n")
  cat("### `Predictive.normal` (normal_inverse_gamma.R)\n")
  cat("- For each data point `x[i]`, it first computes what the posterior parameters (`mu_n`, `kappa_n`, `alpha_n`, `beta_n`) would be if `x[i]` were observed (by calling `PosteriorParameters(mdObj, x[i])`).\n")
  cat("- Then, it uses the formula: `(gamma(alpha_n)/gamma(alpha0)) * ((beta0^alpha0)/(beta_n^alpha_n)) * sqrt(kappa0/kappa_n)`.\n")
  cat("  This corresponds to the density of a non-standardized Student-t distribution, which is the predictive distribution for a Normal likelihood with Normal-Inverse Gamma prior.\n\n")

  cat("### `Predictive.exponential` (exponential_gamma.R)\n")
  cat("- For each `x[i]`, it calculates `alphaPost = alpha0 + 1` (assuming length of x[i] is 1) and `betaPost = beta0 + x[i]`.\n")
  cat("- The predictive likelihood is `(gamma(alphaPost)/gamma(alpha0)) * ((beta0^alpha0)/(betaPost^alphaPost))` which is the density of a Beta-prime (or Gamma-Gamma) distribution scaled appropriately.\n\n")

  cat("### `Predictive.mvnormal` (mvnormal_normal_wishart.R)\n")
  cat("- Calculates posterior parameters `post_params` for each `x[i,]`.\n")
  cat("- The formula involves `pi^(-d/2)`, ratios of `kappa0/kappa_n`, determinants of `Lambda` (prior scale matrix) and `t_n` (posterior scale matrix), and a product of Gamma function ratios. This corresponds to the density of a multivariate Student-t distribution.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `mdObj`: The mixing distribution object, containing `priorParameters` and defining the kernel and base measure.\n")
  cat("- `x`: The data point(s) for which the predictive likelihood is being calculated (numeric vector or matrix).\n")
  cat("- Output: A numeric vector `predictiveArray` of the same length as `x`, containing the predictive likelihoods.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Loop over Data Points**: The functions typically loop through each data point in `x` if multiple are provided.\n")
  cat("2.  **Complexity of Analytical Formula**: The cost per data point depends on the complexity of the analytical formula for the predictive likelihood. This can involve Gamma functions, determinants (for MVN), square roots, and other arithmetic operations.\n")
  cat("3.  **`PosteriorParameters()` Call**: Some implementations (like `Predictive.normal`) call `PosteriorParameters()` internally for each data point, which itself involves calculations based on that point.\n\n")

  cat("## Dependencies\n\n")
  cat("- Base R math functions (`gamma`, `sqrt`, `log`, `det`, `pi`).\n")
  cat("- `PosteriorParameters()` method for the specific conjugate kernel (used by some implementations).\n\n")

  sink()
}

analyze_cluster_label_change <- function() {
  sink("analysis/code_documentation/cluster_label_change.md")

  cat("# ClusterLabelChange() Function Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("The `ClusterLabelChange()` S3 generic method is a utility function called within `ClusterComponentUpdate()`. After a new cluster label (`newLabel`) has been sampled for a data point `i` (which was previously in `currentLabel`), this function updates the state of the Dirichlet Process object (`dpObj`) to reflect this change. This includes updating cluster assignments, counts of points per cluster, the cluster parameters list, and the total number of clusters.\n\n")

  cat("## Implementation Structure\n\n")
  cat("S3 dispatch is used based on the `dpObj`'s conjugacy status:\n")
  cat("- `ClusterLabelChange.conjugate()`\n")
  cat("- `ClusterLabelChange.nonconjugate()`\n")
  cat("The two methods handle slightly different logic regarding how new cluster parameters are incorporated, especially when a new cluster is formed from an auxiliary parameter in the non-conjugate case.\n\n")

  cat("## Algorithmic Steps (Common Logic and `ClusterLabelChange.conjugate`)\n\n")
  cat("Inputs: `dpObj`, data point index `i`, `newLabel`, `currentLabel`, `aux` (auxiliary parameters, primarily for non-conjugate).\n")
  cat("1. Retrieve current state: `pointsPerCluster`, `clusterLabels`, `clusterParameters`, `numLabels`, `mdObj` from `dpObj`.\n")
  cat("2. **Case 1: `newLabel` corresponds to an existing cluster (`newLabel <= numLabels`)**\n")
  cat("   a. Increment count for the new cluster: `pointsPerCluster[newLabel] = pointsPerCluster[newLabel] + 1`.\n")
  cat("   b. Update the data point's label: `clusterLabels[i] = newLabel`.\n")
  cat("   c. Check if the `currentLabel` cluster has become empty (`pointsPerCluster[currentLabel] == 0`):\n")
  cat("      i. If empty: Decrement `numLabels`. Remove the `currentLabel` entry from `pointsPerCluster` and from each parameter array in `clusterParameters`.\n")
  cat("      ii. Re-index `clusterLabels`: Any label greater than `currentLabel` is decremented by 1 to fill the gap.\n")
  cat("3. **Case 2: `newLabel` corresponds to a new cluster being formed (`newLabel > numLabels`)**\n")
  cat("   a. **If `currentLabel` became empty (`pointsPerCluster[currentLabel] == 0`)** (This means the point was the last one in its old cluster, and it's forming a new cluster by itself):\n")
  cat("      i.  **Conjugate**: Draw new parameters for this now re-purposed cluster slot: `post_draw = PosteriorDraw(mdObj, dpObj$data[i, , drop = FALSE])`. Update `clusterParameters[[...]][,,currentLabel]` with `post_draw`.\n")
  cat("      ii. **Non-conjugate**: Assign parameters from the chosen auxiliary parameter: `clusterParameters[[...]][,,currentLabel] = aux[[...]][,, newLabel - numLabels]` (where `newLabel - numLabels` is the index into the `aux` array).\n")
  cat("      iii.Increment its count: `pointsPerCluster[currentLabel] = pointsPerCluster[currentLabel] + 1` (it becomes 1).\n")
  cat("      iv. The data point `i` is already implicitly assigned to `currentLabel` if this path is taken (its label doesn't need to change if it's re-using the slot of its just-emptied cluster for the new parameters).\n")
  cat("   b. **Else (`currentLabel` did not become empty OR it was already a new cluster)**:\n")
  cat("      i.  Assign `clusterLabels[i] = numLabels + 1` (or a re-indexed version if an empty slot was used).\n")
  cat("      ii. Increment `numLabels`.\n")
  cat("      iii.Append 1 to `pointsPerCluster`.\n")
  cat("      iv. **Conjugate**: Draw parameters for the new cluster: `post_draw = PosteriorDraw(mdObj, dpObj$data[i, , drop = FALSE])`. Append `post_draw` to `clusterParameters`.\n")
  cat("      v.  **Non-conjugate**: Take parameters from the chosen auxiliary parameter `aux[[...]][,, newLabel - numLabels]` and append to `clusterParameters`.\n")
  cat("4. Update `dpObj` with the modified `pointsPerCluster`, `clusterLabels`, `clusterParameters`, and `numLabels`.\n\n")

  cat("## `ClusterLabelChange.nonconjugate` Specifics\n\n")
  cat("- The main difference lies in step 3.b.iv/v and 3.a.ii. When a new cluster is created (either by splitting off or by a point moving to an auxiliary parameter slot), its parameters are taken from the `aux` list (the auxiliary parameters drawn in `ClusterComponentUpdate.nonconjugate`). The index `newLabel - numLabels` is used to select which set of auxiliary parameters from `aux` to use.\n\n")

  cat("## Key Data Structures Manipulated\n\n")
  cat("- `dpObj$pointsPerCluster`: Counts are incremented/decremented. May be resized if a cluster is removed or added.\n")
  cat("- `dpObj$clusterLabels`: Assignment for data point `i` is changed.\n")
  cat("- `dpObj$clusterParameters`: Parameters for an existing cluster might be removed (if it becomes empty and its slot is reused or removed). New parameters are added if a new cluster is genuinely created (not reusing an empty slot).\n")
  cat("- `dpObj$numberClusters`: Incremented or decremented.\n")
  cat("- `aux` (Input, non-conjugate only): List of auxiliary parameter sets.\n\n")

  cat("## Performance Considerations\n\n")
  cat("1.  **Array/List Manipulation**: The primary costs come from potentially resizing/modifying `clusterParameters` (a list of arrays) and `pointsPerCluster` (a vector). Appending or removing elements from lists/vectors in R can involve copying data.\n")
  cat("2.  **Re-indexing Labels**: If a cluster is removed, `clusterLabels` needs to be re-indexed, which is a vector operation over `n` elements.\n")
  cat("3.  **`PosteriorDraw` Call (Conjugate, new cluster from empty old one)**: In the conjugate case, if a point leaves its cluster making it empty, and then forms a *new* cluster *by itself reusing that slot*, `PosteriorDraw` is called for that single data point. This is generally fast.\n")
  cat("4.  **Parameter Assignment (Non-conjugate)**: Simply assigns parameters from the `aux` list, which is efficient.\n")
  cat("5.  **Frequency of Calls**: Called once for every data point in every MCMC iteration within `ClusterComponentUpdate`.\n\n")

  cat("## Dependencies\n\n")
  cat("- `PosteriorDraw()` (conjugate case, for initializing a new cluster that reuses an emptied slot, or for initializing a completely new cluster slot).\n")
  cat("- (Implicitly) `PriorDraw()` for generating `aux` parameters in the calling function (`ClusterComponentUpdate.nonconjugate`).\n\n")

  sink()
}

analyze_hierarchical_updates <- function() {
  sink("analysis/code_documentation/hierarchical_updates.md")

  cat("# Hierarchical Dirichlet Process (HDP) Specific Updates Analysis\n\n")

  cat("## Algorithm Overview\n\n")
  cat("For Hierarchical Dirichlet Process models, in addition to updating parameters within each individual DP (`indDP`), global parameters that link these DPs must also be updated. This includes the global base measure $G_0$ (represented by its atoms and stick-breaking weights), the parameters of these global atoms, and the global concentration parameter $\\gamma$.\n\n") # Escaped LaTeX

  cat("Key functions involved:\n")
  cat("- `GlobalParameterUpdate()`: Updates the parameters $\\phi_k$ of the atoms in the global base measure $G_0$.\n") # Escaped LaTeX
  cat("- `UpdateG0()`: Updates the discrete global base measure $G_0$ itself by resampling its stick-breaking weights and potentially atoms. This also updates the mixing distributions of individual DPs (`indDP`) to draw from this new $G_0$.\n")
  cat("- `UpdateGamma()` (via `update_concentration`): Updates the top-level concentration parameter $\\gamma$ for $G_0$.\n\n") # Escaped LaTeX

  cat("## `GlobalParameterUpdate.hierarchical()` (global_parameter_update.R)\n\n")
  cat("### Algorithm Overview\n")
  cat("This function updates the parameters $\\phi_k$ of the shared atoms in the global base measure $G_0$. Each $\\phi_k$ is updated based on all data points from all individual DPs that are currently associated with that specific global atom.\n\n") # Escaped LaTeX

  cat("### Algorithmic Steps\n")
  cat("1. Identify unique global atoms currently in use: `global_labels = unique(unlist(lapply(indDPs, function(dp) match(dp$clusterParameters, G0_atoms))))` (conceptual representation).\n")
  cat("2. For each unique global atom $\\phi_k$ (indexed by `global_labels[i]`):\n") # Escaped LaTeX
  cat("   a. Collect all data points `total_pts` from all individual DPs that are currently assigned to local clusters whose parameters are instances of this global atom $\\phi_k$.\n") # Escaped LaTeX
  cat("      - This involves iterating through each `indDP`.\n")
  cat("      - Identifying which local clusters in `indDP` correspond to the current `global_labels[i]`.\n")
  cat("      - Subsetting the data from `indDP` belonging to these local clusters.\n")
  cat("   b. Draw a new parameter value for this global atom from its posterior distribution given `total_pts` and the hyperprior $H$ (the prior for $G_0$'s atoms): `new_param = PosteriorDraw(mixingDistribution_of_G0, total_pts, ...)`.\n")
  cat("      - The `mixingDistribution` of $G_0$ defines the kernel and prior for $\\phi_k$.\n") # Escaped LaTeX
  cat("      - This step uses `PosteriorDraw`, which will be conjugate or non-conjugate based on $H$.\n")
  cat("   c. Update the `globalParameters` list (which stores the $\\phi_k$ values) with `new_param`.\n") # Escaped LaTeX
  cat("   d. **Crucially**: Update the parameters of all local clusters across all `indDPs` that were instances of the *old* $\\phi_k$ to now be instances of the *new* (updated) $\\phi_k$.\n") # Escaped LaTeX
  cat("3. Update the `mixingDistribution$theta_k` within each `indDP` to reflect the (potentially changed) `globalParameters`.\n\n")

  cat("## `UpdateG0()` (update_g0.R)\n\n")
  cat("### Algorithm Overview\n")
  cat("This function updates the discrete global base measure $G_0 = \\sum \\beta_k \\delta_{\\phi_k}$ by resampling the global stick-breaking weights $\\beta_k$ and using the current set of global atoms $\\phi_k$ (which might have been updated by `GlobalParameterUpdate`). It then updates the individual DPs to draw their cluster parameters from this new $G_0$.\n\n") # Escaped LaTeX

  cat("### Algorithmic Steps\n")
  cat("1. Get current `globalParameters` ($\\{\\phi_k\\}$) and global concentration `gamma`.\n") # Escaped LaTeX
  cat("2. Determine usage of global atoms: Count how many distinct local cluster types (tables) across all `indDPs` are associated with each global atom $\\phi_k$. This forms `globalParamTable$Freq`.\n") # Escaped LaTeX
  cat("3. Resample global stick-breaking weights $\\beta_k$:\n") # Escaped LaTeX
  cat("   a. Use a Dirichlet distribution draw: `dirichlet_draws = gtools::rdirichlet(1, c(globalParamTable$Freq, dpobjlist$gamma))`.\n")
  cat("      The counts `globalParamTable$Freq` are the counts of tables (unique cluster parameter values across all groups) that map to each global atom. `dpobjlist$gamma` is the concentration parameter for the top-level DP generating $G_0$.\n")
  cat("   b. Truncate the stick-breaking process for the 'new' part of $G_0$: Generate `numBreaks` new sticks using `StickBreaking(dpobjlist$gamma + numTables, numBreaks)`, where `numTables` is sum of `globalParamTable$Freq`.\n")
  cat("   c. Combine: `sticks = c(dirichlet_draws[1:numGlobalAtomsUsed], dirichlet_draws[numGlobalAtomsUsed+1] * new_sticks_from_StickBreaking)`.\n")
  cat("4. Generate new atoms if needed for the newly created sticks from the truncation:\n")
  cat("   `priorDraws = PriorDraw(mixingDistribution_of_G0, numBreaks)`.\n")
  cat("5. Combine existing global atoms with newly drawn ones: `postParams` becomes `c(current_global_atoms, priorDraws)`.\n")
  cat("6. For each individual DP (`indDP[[i]]`):\n")
  cat("   a. Generate new local stick-breaking weights `pi_k^(j)` for this DP using the new global sticks `beta_k` (now called `sticks`) and the `indDP[[i]]$alpha`: `newGJ = draw_gj(indDP[[i]]$alpha, sticks)`.\n") # Escaped LaTeX
  cat("   b. Update `indDP[[i]]$mixingDistribution$pi_k = newGJ`.\n")
  cat("   c. Update `indDP[[i]]$mixingDistribution$theta_k = postParams` (all indDPs now point to the same updated set of global atoms and their parameters).\n")
  cat("7. Update `dpobjlist$globalParameters = postParams` and `dpobjlist$globalStick = sticks`.\n\n")

  cat("## `UpdateGamma()` (called via `UpdateAlpha.hierarchical` -> `update_concentration` for the global gamma)\n\n")
  cat("### Algorithm Overview\n")
  cat("Updates the top-level concentration parameter $\\gamma$ for the DP that generates $G_0$, using the same Escobar & West method as `UpdateAlpha.default`, but applied to the 'tables' and 'customers' at the global level.\n\n") # Escaped LaTeX

  cat("### Algorithmic Steps\n")
  cat("1. Count `numParams`: Number of unique global atoms currently used (i.e., number of distinct $\\phi_k$ values that have at least one local cluster mapped to them).\n") # Escaped LaTeX
  cat("2. Count `numTables`: Total number of distinct local cluster types across all `indDPs` (i.e., sum of the number of unique cluster parameters in each group, before collapsing to global atoms, effectively this is the sum of `globalParamTable$Freq` from `UpdateG0`).\n")
  cat("3. Call `update_concentration(current_gamma, numTables, numParams, gammaPriors)`.\n\n")

  cat("## Key Data Structures Used\n\n")
  cat("- `dpobjlist$indDP`: List of individual Dirichlet Process objects.\n")
  cat("- `dpobjlist$globalParameters`: List of arrays storing the parameters $\\{\\phi_k\\}$ of the global atoms.\n") # Escaped LaTeX
  cat("- `dpobjlist$globalStick`: Vector of stick-breaking weights $\\{\\beta_k\\}$ for $G_0$.\n") # Escaped LaTeX
  cat("- `dpobjlist$gamma`: Global concentration parameter.\n")
  cat("- `dpobjlist$gammaPriors`: Prior for `gamma`.\n")
  cat("- `indDP[[i]]$clusterParameters`: Parameters of local clusters.\n")
  cat("- `indDP[[i]]$clusterLabels`: Local cluster assignments.\n")
  cat("- `indDP[[i]]$mixingDistribution$pi_k`: Local stick-breaking weights (pointing to global atoms).\n")
  cat("- `indDP[[i]]$mixingDistribution$theta_k`: Reference to global atoms $\\{\\phi_k\\}$.\n\n") # Escaped LaTeX

  cat("## Performance Considerations\n\n")
  cat("1.  **`GlobalParameterUpdate`**:\n")
  cat("    -   Iterates over unique global atoms.\n")
  cat("    -   For each global atom, iterates over all `indDPs` and their local clusters to gather data.\n")
  cat("    -   `PosteriorDraw()` for each global atom can be expensive, especially if non-conjugate and data pooled from many groups is large.\n")
  cat("2.  **`UpdateG0`**:\n")
  cat("    -   Counting global atom usage involves iterating through `indDPs` and their cluster parameters.\n")
  cat("    -   `gtools::rdirichlet` and `StickBreaking` are generally efficient.\n")
  cat("    -   `PriorDraw` for new global atoms depends on `numBreaks`.\n")
  cat("    -   The loop to update `pi_k` for each `indDP` involves `draw_gj`, which itself has a loop.\n")
  cat("3.  **Data Aggregation**: Collecting all data points associated with a global atom (`total_pts` in `GlobalParameterUpdate`) can be memory and time-consuming if not done efficiently.\n\n")

  cat("## Dependencies\n\n")
  cat("- All core DP functions for individual DPs (`PosteriorDraw`, `PriorDraw`, `Likelihood` for the kernel of $\\phi_k$).\n") # Escaped LaTeX
  cat("- `gtools::rdirichlet`.\n")
  cat("- `StickBreaking()`, `draw_gj()`.\n")
  cat("- `update_concentration()` (for `UpdateGamma`).\n")
  cat("- `true_cluster_labels()` for correctly matching multi-dimensional parameters.\n\n")

  sink()
}


# --- Call the main documentation function ---
document_algorithms()

cat("Algorithm analysis markdown files generated in analysis/code_documentation/\n")
