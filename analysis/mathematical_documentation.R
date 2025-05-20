cat("DEBUG: Defining function document_all_maths...\n")
# Main function to call all documentation generators
document_all_maths <- function() {
  # Print a message to confirm script execution start
  cat("Starting mathematical_documentation_generator.R script execution...\n")

  base_dir <- "analysis/code_documentation" # Define base directory for outputs
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
    cat("Created directory:", base_dir, "\n")
  }

  document_mathematical_basis(file.path(base_dir, "mathematical_basis.md"))
  document_likelihood_math(base_dir)
  document_prior_draw_math(base_dir)
  document_posterior_draw_math(base_dir)
  document_predictive_likelihood_math(base_dir)
  document_mcmc_components_math(base_dir) # For CRP, MH details
  document_concentration_param_math(file.path(base_dir, "concentration_parameter_update.md"))
  document_hierarchical_math(base_dir) # For HDP specific updates

  cat("All mathematical documentation generation complete.\nFiles are in", base_dir, "\n")
}
cat("DEBUG: Defining function document_mathematical_basis...\n")
# 1. General Mathematical Basis (from user input, with corrections)
document_mathematical_basis <- function(filepath) {
  cat("DEBUG: Executing document_mathematical_basis for", filepath, "\n")
  sink(filepath)

  writeLines("# Mathematical Foundation of dirichletprocess Algorithms

## Dirichlet Process Mixture Model (DPMM)

The model assumes data $y_i$ are generated from a mixture model where the mixing distribution $G$ itself is drawn from a Dirichlet Process (DP). The hierarchical representation is:

$$y_i | \\theta_i \\sim F(\\cdot | \\theta_i)$$
$$\\theta_i | G \\sim G$$
$$G | G_0, \\alpha \\sim DP(G_0, \\alpha)$$

Where:
- $y_i$ is the $i$-th observation.
- $F(\\cdot | \\theta_i)$ is the likelihood function for an observation given parameters $\\theta_i$ (e.g., Gaussian $N(y_i | \\mu_i, \\sigma_i^2)$, Beta $Beta(y_i | a_i, b_i)$).
- $\\theta_i$ are the parameters specific to the $i$-th observation.
- $G$ is a random probability measure (the mixing distribution) drawn from a Dirichlet Process.
- $G_0$ is the base measure, representing the prior belief about the distribution of parameters $\\theta$. It is the expected value of $G$, $E[G] = G_0$.
- $\\alpha > 0$ is the concentration parameter. Larger $\\alpha$ leads to $G$ being closer to $G_0$ and implies a larger number of distinct parameter values (clusters) a priori.

Due to the discreteness of draws from a DP, many $\\theta_i$ will share common values. Let $\\{\\phi_k\\}_{k=1}^K$ be the set of unique parameter values, where $K$ is the number of clusters.

## Chinese Restaurant Process (CRP) / Gibbs Sampling for Cluster Assignments (Neal's Algorithm 3)

The Chinese Restaurant Process (CRP) formulation provides a constructive way to understand the clustering property of the DP. For Gibbs sampling, we integrate out $G$. The conditional posterior probability for assigning data point $y_i$ to cluster $k$ (parameterized by $\\phi_k$), given all other assignments $c_{-i}$ and current cluster parameters, is (Neal, 2000, Algorithm 3):

$$P(c_i = k | c_{-i}, \\mathbf{y}, \\alpha, G_0, \\{\\phi_j\\}_{j \\neq k}) \\propto \\begin{cases}
n_{k,-i} \\cdot F(y_i | \\phi_k) & \\text{for an existing cluster } k \\text{ (where } n_{k,-i} > 0) \\\\
\\alpha \\cdot \\int F(y_i | \\phi) dG_0(\\phi) & \\text{for a new cluster (i.e., } k = K_{new})
\\end{cases}$$

Where:
- $c_i$ is the cluster assignment for data point $y_i$.
- $c_{-i}$ denotes the cluster assignments for all data points except $y_i$.
- $n_{k,-i}$ is the number of data points in cluster $k$, excluding point $y_i$.
- $F(y_i | \\phi_k)$ is the likelihood of data point $y_i$ given the parameters $\\phi_k$ of cluster $k$.
- $\\alpha$ is the concentration parameter.
- $\\int F(y_i | \\phi) dG_0(\\phi)$ is the marginal likelihood of $y_i$ under the base measure $G_0$, often called the predictive likelihood. This is used when proposing a new cluster, whose parameters $\\phi_{K_{new}}$ would be drawn from $G_0$ (or a posterior based on $y_i$ if conjugate).

This is Neal's (2000) Algorithm 3. Algorithm 2 is similar but marginalizes over $\\phi_k$ as well, which is common for fully conjugate models.

## Cluster Parameter Updates (Neal's Algorithm 3)

After re-assigning all $c_i$, the parameters $\\phi_k$ for each cluster $k$ (now containing points $D_k = \\{y_i : c_i = k\\}$) are updated.

### Conjugate Case

If $G_0$ is conjugate to the likelihood $F$, the parameters $\\phi_k$ for each cluster $k$ are drawn from their posterior distribution:

$$\\phi_k | D_k, G_0 \\sim p(\\phi_k | D_k, G_0)$$

This posterior is analytically known. For example, if $F$ is Gaussian and $G_0$ is Normal-Inverse-Gamma, the posterior $p(\\phi_k | D_k, G_0)$ is also Normal-Inverse-Gamma with updated hyperparameters.

### Non-conjugate Case (Neal's Algorithm 8 for assignments, Metropolis-Hastings for parameters)

When $G_0$ is not conjugate to $F$, sampling $\\phi_k$ directly from the posterior is not possible. Instead, MCMC methods like Metropolis-Hastings are used. The target density is:

$$p(\\phi_k | D_k, G_0) \\propto \\left( \\prod_{y_j \\in D_k} F(y_j | \\phi_k) \\right) p_0(\\phi_k)$$

where $p_0(\\phi_k)$ is the prior density of $\\phi_k$ from $G_0$.
For cluster assignments in non-conjugate cases, Neal's Algorithm 8 is often used, which introduces $m$ auxiliary parameters drawn from $G_0$ to help propose new clusters, avoiding the direct calculation of $\\int F(y_i | \\phi) dG_0(\\phi)$.

## Concentration Parameter $\\alpha$ Update (Escobar & West, 1995)

The concentration parameter $\\alpha$ can be given a prior (e.g., Gamma distribution $p(\\alpha) = \\text{Gamma}(\\alpha | a_\\alpha, b_\\alpha)$) and updated using an auxiliary variable method (Escobar and West, 1995). Given $K$ current clusters and $N$ data points, the update involves drawing an auxiliary variable $\\eta \\sim \\text{Beta}(\\alpha+1, N)$. The posterior for $\\alpha$ is then a mixture of two Gamma distributions:

$$p(\\alpha | K, N, \\eta, \\text{priors}) = \\pi_\\eta \\text{Gamma}(\\alpha | a_\\alpha+K, b_\\alpha-\\log(\\eta)) + (1-\\pi_\\eta) \\text{Gamma}(\\alpha | a_\\alpha+K-1, b_\\alpha-\\log(\\eta))$$

Where the mixture weight $\\pi_\\eta$ depends on $a_\\alpha, K, N, b_\\alpha,$ and $\\log(\\eta)$.
")

  sink()
  cat("Generated:", filepath, "\n")
}
cat("DEBUG: Defining function document_likelihood_math...\n")
# 2. Likelihood Functions
document_likelihood_math <- function(base_dir) {
  cat("DEBUG: Executing document_likelihood_math...\n")
  # Gaussian Likelihood
  filepath_gaussian <- file.path(base_dir, "likelihood_gaussian.md")
  cat("DEBUG: Writing to (Gaussian):", filepath_gaussian, "\n")
  sink(filepath_gaussian)

  writeLines("# Gaussian Likelihood Function $F(y | \\mu, \\sigma^2)$

## Mathematical Formula

The likelihood for a data point $y$ given mean $\\mu$ and variance $\\sigma^2$ (or standard deviation $\\sigma$) is:

$$F(y | \\mu, \\sigma^2) = \\frac{1}{\\sqrt{2\\pi\\sigma^2}} \\exp\\left(-\\frac{(y - \\mu)^2}{2\\sigma^2}\\right) = \\frac{1}{\\sigma\\sqrt{2\\pi}} \\exp\\left(-\\frac{(y - \\mu)^2}{2\\sigma^2}\\right)$$

## Implementation Considerations

1. For numerical stability, compute log-likelihood first:
   $$\\log F(y | \\mu, \\sigma^2) = -\\frac{1}{2}\\log(2\\pi\\sigma^2) - \\frac{(y - \\mu)^2}{2\\sigma^2} = -\\log(\\sigma) - \\frac{1}{2}\\log(2\\pi) - \\frac{(y - \\mu)^2}{2\\sigma^2}$$

2. Vectorization: Compute for multiple data points $y_i$ simultaneously given a single $(\\mu, \\sigma^2)$ or for a single $y_i$ against multiple $(\\mu_k, \\sigma_k^2)$.

3. Ensure $\\sigma^2 > 0$ (or $\\sigma > 0$). Handle cases where $\\sigma$ is very small to avoid overflow/underflow.")

  sink()
  cat("Generated:", filepath_gaussian, "\n")

  # Beta Likelihood
  filepath_beta <- file.path(base_dir, "likelihood_beta.md")
  cat("DEBUG: Writing to (Beta):", filepath_beta, "\n")
  sink(filepath_beta)

  writeLines("# Beta Likelihood Function $F(y | \\mu, \\tau, M_{max})$

## Mathematical Formula

The `dirichletprocess` package often parameterizes the Beta distribution for data $y \\in [0, M_{max}]$ using a mean $\\mu_y \\in [0, M_{max}]$ and a precision parameter $\\tau > 0$. Let $x = y/M_{max}$ be the scaled data point in $[0,1]$, with mean $\\mu_x = \\mu_y/M_{max}$. The standard Beta distribution for $x$ is $\\text{Beta}(x | a, b)$. The parameters $a$ and $b$ are related to $\\mu_x$ and $\\tau$ by:

$$a = \\mu_x \\tau$$
$$b = (1 - \\mu_x) \\tau$$

The likelihood for $y$ is then given by the transformed density:

$$F(y | \\mu_y, \\tau, M_{max}) = \\frac{1}{M_{max}} \\text{Beta}\\left(\\frac{y}{M_{max}} \\Big| a = \\frac{\\mu_y}{M_{max}}\\tau, b = \\left(1 - \\frac{\\mu_y}{M_{max}}\\right)\\tau \\right)$$
$$= \\frac{1}{M_{max}} \\frac{\\Gamma(a+b)}{\\Gamma(a)\\Gamma(b)} \\left(\\frac{y}{M_{max}}\\right)^{a-1} \\left(1 - \\frac{y}{M_{max}}\\right)^{b-1}$$

Where $0 \\le y \\le M_{max}$, $0 < \\mu_y < M_{max}$, and $\\tau > 0$.

## Implementation Considerations

1. Ensure $a > 0$ and $b > 0$. This implies $\\mu_y \\in (0, M_{max})$ and $\\tau > 0$.
2. Handle edge cases: If $y=0$ or $y=M_{max}$, the likelihood can be 0 or infinity if $a<1$ or $b<1$. Standard library functions for Beta PDF usually handle this.
3. Log-likelihood is preferred for products: $\\log F(y | \\dots) = -\\log(M_{max}) + \\log\\Gamma(a+b) - \\log\\Gamma(a) - \\log\\Gamma(b) + (a-1)\\log(y/M_{max}) + (b-1)\\log(1-y/M_{max})$.")

  sink()
  cat("Generated:", filepath_beta, "\n")

  # Weibull Likelihood
  filepath_weibull <- file.path(base_dir, "likelihood_weibull.md")
  cat("DEBUG: Writing to (Weibull):", filepath_weibull, "\n")
  sink(filepath_weibull)

  writeLines("# Weibull Likelihood Function $F(y | k, \\lambda)$

## Mathematical Formula

The Weibull distribution is often parameterized by a shape parameter $k > 0$ and a scale parameter $\\lambda > 0$. The probability density function is:

$$F(y | k, \\lambda) = \\frac{k}{\\lambda} \\left(\\frac{y}{\\lambda}\\right)^{k-1} \\exp\\left(-\\left(\\frac{y}{\\lambda}\\right)^k\\right)$$

for $y \\ge 0$.

The R `dweibull(x, shape, scale)` function uses this parameterization directly.

## Implementation Considerations

1. Ensure $y \\ge 0$, $k > 0$, and $\\lambda > 0$.
2. Log-likelihood:
   $$\\log F(y | k, \\lambda) = \\log(k) - \\log(\\lambda) + (k-1)(\\log(y) - \\log(\\lambda)) - \\left(\\frac{y}{\\lambda}\\right)^k$$
3. Be careful with $y=0$: If $k=1$ (Exponential), $f(0|1,\\lambda)=1/\\lambda$. If $k>1$, $f(0|k,\\lambda)=0$. If $0<k<1$, $f(0|k,\\lambda)=\\infty$.")

  sink()
  cat("Generated:", filepath_weibull, "\n")

  # Exponential Likelihood
  filepath_exp <- file.path(base_dir, "likelihood_exponential.md")
  cat("DEBUG: Writing to (Exponential):", filepath_exp, "\n")
  sink(filepath_exp)

  writeLines("# Exponential Likelihood Function $F(y | \\lambda_{rate})$

## Mathematical Formula

The Exponential distribution is parameterized by a rate parameter $\\lambda_{rate} > 0$. The probability density function is:

$$F(y | \\lambda_{rate}) = \\lambda_{rate} e^{-\\lambda_{rate} y}$$

for $y \\ge 0$.

## Implementation Considerations

1. Ensure $y \\ge 0$ and $\\lambda_{rate} > 0$.
2. Log-likelihood: $\\log F(y | \\lambda_{rate}) = \\log(\\lambda_{rate}) - \\lambda_{rate} y$.")

  sink()
  cat("Generated:", filepath_exp, "\n")

  # Multivariate Normal Likelihood
  filepath_mvn <- file.path(base_dir, "likelihood_mvnormal.md")
  cat("DEBUG: Writing to (MVN):", filepath_mvn, "\n")
  sink(filepath_mvn)

  writeLines("# Multivariate Normal Likelihood Function $F(\\mathbf{y} | \\boldsymbol{\\mu}, \\boldsymbol{\\Sigma})$

## Mathematical Formula

For a $d$-dimensional data point $\\mathbf{y}$, given a mean vector $\\boldsymbol{\\mu}$ ($d \\times 1$) and a covariance matrix $\\boldsymbol{\\Sigma}$ ($d \\times d$, symmetric positive definite), the likelihood is:

$$F(\\mathbf{y} | \\boldsymbol{\\mu}, \\boldsymbol{\\Sigma}) = (2\\pi)^{-d/2} |\\boldsymbol{\\Sigma}|^{-1/2} \\exp\\left(-\\frac{1}{2}(\\mathbf{y} - \\boldsymbol{\\mu})^T \\boldsymbol{\\Sigma}^{-1} (\\mathbf{y} - \\boldsymbol{\\mu})\\right)$$

Where $|\\cdot|$ denotes the determinant and $T$ denotes transpose.

## Implementation Considerations

1. Log-likelihood is numerically preferred:
   $$\\log F(\\mathbf{y} | \\boldsymbol{\\mu}, \\boldsymbol{\\Sigma}) = -\\frac{d}{2}\\log(2\\pi) - \\frac{1}{2}\\log|\\boldsymbol{\\Sigma}| - \\frac{1}{2}(\\mathbf{y} - \\boldsymbol{\\mu})^T \\boldsymbol{\\Sigma}^{-1} (\\mathbf{y} - \\boldsymbol{\\mu})$$

2. The term $(\\mathbf{y} - \\boldsymbol{\\mu})^T \\boldsymbol{\\Sigma}^{-1} (\\mathbf{y} - \\boldsymbol{\\mu})$ is the Mahalanobis distance squared.
3. Efficient computation of $\\log|\\boldsymbol{\\Sigma}|$ and $\\boldsymbol{\\Sigma}^{-1}$ (or solving $\\boldsymbol{\\Sigma}\\mathbf{x} = (\\mathbf{y} - \\boldsymbol{\\mu})$ for $\\mathbf{x}$ and then computing $\\mathbf{x}^T(\\mathbf{y} - \\boldsymbol{\\mu})$) is crucial. Often, the Cholesky decomposition of $\\boldsymbol{\\Sigma} = LL^T$ is used:
   - $\\log|\\boldsymbol{\\Sigma}| = 2 \\sum_j \\log(L_{jj})$
   - The quadratic form can be computed by solving $L\\mathbf{z} = (\\mathbf{y} - \\boldsymbol{\\mu})$ for $\\mathbf{z}$, then the quadratic form is $\\mathbf{z}^T\\mathbf{z}$.
4. Ensure $\\boldsymbol{\\Sigma}$ is symmetric and positive definite.")

  sink()
  cat("Generated:", filepath_mvn, "\n")
}
cat("DEBUG: Defining function document_prior_draw_math...\n")
# 3. Prior Draw Functions
document_prior_draw_math <- function(base_dir) {
  cat("DEBUG: Executing document_prior_draw_math...\n")
  # Prior Draw for Gaussian (Normal-Inverse Gamma)
  filepath_norm_ig <- file.path(base_dir, "prior_draw_gaussian_norm_invgamma.md")
  cat("DEBUG: Writing to", filepath_norm_ig, "\n")
  sink(filepath_norm_ig)

  writeLines("# Prior Draw: Gaussian Parameters from Normal-Inverse Gamma ($G_0$)

For a Gaussian likelihood $N(y | \\mu, \\sigma^2)$, a common conjugate base measure $G_0$ is the Normal-Inverse Gamma distribution for $(\\mu, \\sigma^2)$.

Hyperparameters for $G_0$: $\\mu_0, \\kappa_0 > 0, a_0 > 0, b_0 > 0$.

Sampling procedure:
1. Draw variance $\\sigma^2$ from an Inverse Gamma distribution:
   $$\\sigma^2 \\sim \\text{Inv-Gamma}(a_0, b_0)$$
   (This is equivalent to drawing precision $\\lambda = 1/\\sigma^2 \\sim \\text{Gamma}(a_0, b_0)$).
2. Draw mean $\\mu$ conditional on the drawn $\\sigma^2$ from a Normal distribution:
   $$\\mu | \\sigma^2 \\sim N(\\mu_0, \\sigma^2 / \\kappa_0)$$

The pair $(\\mu, \\sigma^2)$ is one draw from $G_0$.")

  sink()
  cat("Generated:", filepath_norm_ig, "\n")

  # Prior Draw for Beta (e.g., Uniform for mean, Gamma for precision/scale)
  filepath_beta_prior <- file.path(base_dir, "prior_draw_beta.md")
  cat("DEBUG: Writing to", filepath_beta_prior, "\n")
  sink(filepath_beta_prior)

  writeLines("# Prior Draw: Beta Parameters from $G_0$

For a Beta likelihood $F(y | \\mu_y, \\tau, M_{max})$ (mean $\\mu_y \\in [0, M_{max}]$, precision $\\tau$), a typical non-conjugate $G_0$ might define independent priors for $\\mu_y$ and $\\tau$.

Hyperparameters for $G_0$ could be:
- For $\\mu_y$: Often $U(0, M_{max})$.
- For $\\tau$: e.g., $\\text{Gamma}(a_\\tau, b_\\tau)$ or $\\text{Inv-Gamma}(a_\\tau, b_\\tau)$. The `dirichletprocess` R code for `beta_uniform_gamma` uses an Inverse Gamma prior for `nu` (scale parameter, $\\tau = \\text{nu}$ in that context), meaning `nu` ~ InvGamma(alpha0, beta0).

Sampling procedure (example using $U(0,M_{max})$ for $\\mu_y$ and $\\text{Gamma}(a_\\tau,b_\\tau)$ for $\\tau$ as precision):
1. Draw mean $\\mu_y$:
   $$\\mu_y \\sim U(0, M_{max})$$
2. Draw precision $\\tau$:
   $$\\tau \\sim \\text{Gamma}(a_\\tau, b_\\tau)$$

The pair $(\\mu_y, \\tau)$ is one draw from $G_0$.")

  sink()
  cat("Generated:", filepath_beta_prior, "\n")

  # Prior Draw for Weibull
  filepath_weibull_prior <- file.path(base_dir, "prior_draw_weibull.md")
  cat("DEBUG: Writing to", filepath_weibull_prior, "\n")
  sink(filepath_weibull_prior)

  writeLines("# Prior Draw: Weibull Parameters ($k, \\lambda$) from $G_0$

For a Weibull likelihood $F(y | k, \\lambda)$ (shape $k$, scale $\\lambda$), $G_0$ might define independent priors for $k$ and $\\lambda$.

Hyperparameters for $G_0$ (example from `weibull_uniform_gamma.R`):
- For shape $k$: $U(0, \\phi_{max})$ where $\\phi_{max}$ is `priorParameters[1]`.
- For scale $\\lambda$: $\\text{Inv-Gamma}(a_\\lambda, b_\\lambda)$ where $a_\\lambda$=`priorParameters[2]`, $b_\\lambda$=`priorParameters[3]`.

Sampling procedure:
1. Draw shape $k$:
   $$k \\sim U(0, \\phi_{max})$$
2. Draw scale $\\lambda$:
   $$\\lambda \\sim \\text{Inv-Gamma}(a_\\lambda, b_\\lambda)$$
   (Equivalent to $1/\\lambda \\sim \\text{Gamma}(a_\\lambda, b_\\lambda)$).

The pair $(k, \\lambda)$ is one draw from $G_0$.")

  sink()
  cat("Generated:", filepath_weibull_prior, "\n")

  # Prior Draw for Exponential
  filepath_exp_prior <- file.path(base_dir, "prior_draw_exponential.md")
  cat("DEBUG: Writing to", filepath_exp_prior, "\n")
  sink(filepath_exp_prior)

  writeLines("# Prior Draw: Exponential Parameter (rate $\\lambda_{rate}$) from $G_0$

For an Exponential likelihood $F(y | \\lambda_{rate})$, a common conjugate base measure $G_0$ for the rate parameter $\\lambda_{rate}$ is a Gamma distribution.

Hyperparameters for $G_0$: shape $a_0 > 0$, rate $b_0 > 0$.

Sampling procedure:
1. Draw rate $\\lambda_{rate}$ from a Gamma distribution:
   $$\\lambda_{rate} \\sim \\text{Gamma}(a_0, b_0)$$

This $\\lambda_{rate}$ is one draw from $G_0$.")

  sink()
  cat("Generated:", filepath_exp_prior, "\n")

  # Prior Draw for MVN (Normal-Wishart)
  filepath_mvn_prior <- file.path(base_dir, "prior_draw_mvnormal_norm_wishart.md")
  cat("DEBUG: Writing to", filepath_mvn_prior, "\n")
  sink(filepath_mvn_prior)

  writeLines("# Prior Draw: Multivariate Normal Parameters from Normal-Wishart ($G_0$)

For an MVN likelihood $N_d(\\mathbf{y} | \\boldsymbol{\\mu}, \\boldsymbol{\\Sigma})$, a conjugate base measure $G_0$ for $(\\boldsymbol{\\mu}, \\boldsymbol{\\Lambda} = \\boldsymbol{\\Sigma}^{-1})$ is Normal-Wishart.

Hyperparameters for $G_0$: mean vector $\\boldsymbol{m}_0$, scalar $\\kappa_0 > 0$, degrees of freedom $\\nu_0 > d-1$, scale matrix $\\mathbf{S}_0$ ($d \\times d$, symmetric positive definite for $\\boldsymbol{\\Lambda}$).

Sampling procedure:
1. Draw precision matrix $\\boldsymbol{\\Lambda}$ from a Wishart distribution:
   $$\\boldsymbol{\\Lambda} \\sim \\text{Wishart}(\\nu_0, \\mathbf{S}_0)$$
   (Note: R's `rWishart(nu, V)` samples from Wishart with df $\\nu$ and scale matrix $V$. If $\\mathbf{S}_0$ is the scale matrix for $\\boldsymbol{\\Lambda}$, this is direct. Sometimes $\\mathbf{S}_0$ is parameterized as $(E[\\boldsymbol{\\Lambda}])^{-1}/\\nu_0$).

2. Draw mean vector $\\boldsymbol{\\mu}$ conditional on $\\boldsymbol{\\Lambda}$ from a Multivariate Normal distribution:
   $$\\boldsymbol{\\mu} | \\boldsymbol{\\Lambda} \\sim N_d(\\boldsymbol{m}_0, (\\kappa_0 \\boldsymbol{\\Lambda})^{-1})$$

The pair $(\\boldsymbol{\\mu}, \\boldsymbol{\\Sigma} = \\boldsymbol{\\Lambda}^{-1})$ is one draw from $G_0$.")

  sink()
  cat("Generated:", filepath_mvn_prior, "\n")
}
cat("DEBUG: Defining function document_posterior_draw_math...\n")
# 4. Posterior Draw Functions
document_posterior_draw_math <- function(base_dir) {
  cat("DEBUG: Executing document_posterior_draw_math...\n")
  # Posterior Draw for Gaussian (Normal-Inverse Gamma)
  filepath_norm_ig_post <- file.path(base_dir, "posterior_draw_gaussian_norm_invgamma.md")
  cat("DEBUG: Writing to", filepath_norm_ig_post, "\n")
  sink(filepath_norm_ig_post)

  writeLines("# Posterior Draw: Gaussian Parameters (Normal-Inverse Gamma Conjugacy)

Given data $D_k = \\{y_j : c_j=k\\}$ assigned to cluster $k$, with $n_k = |D_k|$, sample mean $\\bar{y}_k = \\frac{1}{n_k}\\sum_{y_j \\in D_k} y_j$, and sum of squared deviations $S_k = \\sum_{y_j \\in D_k} (y_j - \\bar{y}_k)^2$.

Prior hyperparameters (from $G_0$): $\\mu_0, \\kappa_0, a_0, b_0$.

Posterior hyperparameters:
$$\\kappa_n = \\kappa_0 + n_k$$
$$\\mu_n = \\frac{\\kappa_0 \\mu_0 + n_k \\bar{y}_k}{\\kappa_n}$$
$$a_n = a_0 + \\frac{n_k}{2}$$
$$b_n = b_0 + \\frac{1}{2}S_k + \\frac{\\kappa_0 n_k (\\bar{y}_k - \\mu_0)^2}{2(\\kappa_0 + n_k)}$$

Sampling procedure for cluster parameters $(\\mu_k, \\sigma_k^2)$:
1. Draw variance $\\sigma_k^2$ from the posterior Inverse Gamma distribution:
   $$\\sigma_k^2 | D_k, G_0 \\sim \\text{Inv-Gamma}(a_n, b_n)$$
2. Draw mean $\\mu_k$ conditional on $\\sigma_k^2$ from the posterior Normal distribution:
   $$\\mu_k | \\sigma_k^2, D_k, G_0 \\sim N(\\mu_n, \\sigma_k^2 / \\kappa_n)$$")

  sink()
  cat("Generated:", filepath_norm_ig_post, "\n")

  # Posterior Draw for Beta (Non-conjugate)
  filepath_beta_post <- file.path(base_dir, "posterior_draw_beta_nonconjugate.md")
  cat("DEBUG: Writing to", filepath_beta_post, "\n")
  sink(filepath_beta_post)

  writeLines("# Posterior Draw: Beta Parameters (Non-Conjugate via Metropolis-Hastings)

For a Beta likelihood $F(y | \\mu_y, \\tau, M_{max})$ with non-conjugate priors on $(\\mu_y, \\tau)$ (e.g., $U(0,M_{max})$ for $\\mu_y$, $\\text{Gamma}(a_\\tau, b_\\tau)$ for $\\tau$), we use Metropolis-Hastings to sample from the posterior.

Given data $D_k = \\{y_j : c_j=k\\}$ for cluster $k$.
Target posterior density:
$$p(\\mu_y, \\tau | D_k, G_0) \\propto \\left( \\prod_{y_j \\in D_k} F(y_j | \\mu_y, \\tau, M_{max}) \\right) p_0(\\mu_y) p_0(\\tau)$$
where $p_0(\\cdot)$ are the prior densities from $G_0$.

Metropolis-Hastings steps for $(\\mu_y, \\tau)$:
1. Initialize $(\\mu_y^{(t)}, \\tau^{(t)})$.
2. Propose $(\\mu_y', \\tau')$ from $q((\\mu_y', \\tau') | (\\mu_y^{(t)}, \\tau^{(t)}))$. E.g., independent random walks for $\\mu_y$ (clipped to $[0, M_{max}]$) and $\\log(\\tau)$ or $\\tau$ (ensuring positivity).
3. Calculate acceptance ratio $\\mathcal{A}$ (see `metropolis_hastings.md` for general formula).
4. Accept or reject the proposal to get $(\\mu_y^{(t+1)}, \\tau^{(t+1)})$.

The `dirichletprocess` package uses specific proposal mechanisms (`MhParameterProposal.beta`) and prior densities (`PriorDensity.beta`) within the MH algorithm.")

  sink()
  cat("Generated:", filepath_beta_post, "\n")

  # Posterior Draw for Weibull (Non-conjugate)
  filepath_weibull_post <- file.path(base_dir, "posterior_draw_weibull_nonconjugate.md")
  cat("DEBUG: Writing to", filepath_weibull_post, "\n")
  sink(filepath_weibull_post)

  writeLines("# Posterior Draw: Weibull Parameters (Non-Conjugate via Metropolis-Hastings with Gibbs Step)

For a Weibull likelihood $F(y | k, \\lambda)$ (shape $k$, scale $\\lambda$) with potentially non-conjugate priors (e.g., $U(0, \\phi_{max})$ for $k$, $\\text{Inv-Gamma}(a_\\lambda,b_\\lambda)$ for $\\lambda$).

Given data $D_k = \\{y_j : c_j=k\\}$ for cluster $k$.
Target posterior density:
$$p(k, \\lambda | D_k, G_0) \\propto \\left( \\prod_{y_j \\in D_k} F(y_j | k, \\lambda) \\right) p_0(k) p_0(\\lambda)$$

The `dirichletprocess` package often uses a Metropolis-within-Gibbs approach or a specialized MH for Weibull:
1. Sample shape $k'$ from a proposal $q(k'|k^{(t)})$ (e.g., random walk on $k$ or $\\log k$).
2. Given $k'$, sample scale $\\lambda'$ from its full conditional posterior $p(\\lambda | k', D_k, G_0)$. If the prior $p_0(\\lambda)$ is Inv-Gamma (i.e., $1/\\lambda \\sim \\text{Gamma}(a_\\lambda, b_\\lambda)$), then the full conditional for $1/\\lambda$ is:
   $$1/\\lambda | k', D_k, G_0 \\sim \\text{Gamma}\\left(a_\\lambda + n_k, b_\\lambda + \\sum_{y_j \\in D_k} (y_j)^{k'}\\right)$$
   (The R code for `PosteriorDraw.weibull` samples $1/\\lambda$ from $\\text{Gamma}(n_k + a_0, \\sum y_j^{k'} + b_0)$ implying $a_0, b_0$ are prior params for $1/\\lambda$. The $(y_j/M_{scale})^{k'}$ form is more general if a scaling factor $M_{scale}$ were involved, but typically not for Weibull scale parameter $\\lambda$ itself in this context).
3. Calculate acceptance ratio for $k'$ based on the marginal likelihood of $k'$ (integrating out $\\lambda$) or using the joint proposal $(k', \\lambda')$. The R code `MetropolisHastings.weibull` uses the joint proposal where $\\lambda'$ is drawn from its conditional given $k'$.")

  sink()
  cat("Generated:", filepath_weibull_post, "\n")

  # Posterior Draw for Exponential (Gamma Conjugacy)
  filepath_exp_post <- file.path(base_dir, "posterior_draw_exponential_gamma.md")
  cat("DEBUG: Writing to", filepath_exp_post, "\n")
  sink(filepath_exp_post)

  writeLines("# Posterior Draw: Exponential Parameter (Gamma Conjugacy)

For an Exponential likelihood $F(y | \\lambda_{rate})$ and a Gamma prior $G_0 = \\text{Gamma}(\\lambda_{rate} | a_0, b_0)$ on the rate parameter.

Given data $D_k = \\{y_j : c_j=k\\}$ for cluster $k$, with $n_k = |D_k|$ and sum $S_k = \\sum_{y_j \\in D_k} y_j$.

Posterior hyperparameters for $\\lambda_{rate,k}$:
$$a_n = a_0 + n_k$$
$$b_n = b_0 + S_k$$

Sampling procedure for cluster rate parameter $\\lambda_{rate,k}$:
1. Draw $\\lambda_{rate,k}$ from the posterior Gamma distribution:
   $$\\lambda_{rate,k} | D_k, G_0 \\sim \\text{Gamma}(a_n, b_n)$$")

  sink()
  cat("Generated:", filepath_exp_post, "\n")

  # Posterior Draw for MVN (Normal-Wishart)
  filepath_mvn_post <- file.path(base_dir, "posterior_draw_mvnormal_norm_wishart.md")
  cat("DEBUG: Writing to", filepath_mvn_post, "\n")
  sink(filepath_mvn_post)

  writeLines("# Posterior Draw: MVN Parameters (Normal-Wishart Conjugacy)

Given data $D_k = \\{\\mathbf{y}_j : c_j=k\\}$ for cluster $k$, $n_k = |D_k|$, sample mean $\\bar{\\mathbf{y}}_k = \\frac{1}{n_k}\\sum \\mathbf{y}_j$, scatter matrix $\\mathbf{S}_{scat,k} = \\sum (\\mathbf{y}_j - \\bar{\\mathbf{y}}_k)(\\mathbf{y}_j - \\bar{\\mathbf{y}}_k)^T$.

Prior hyperparameters ($G_0$): $\\boldsymbol{m}_0, \\kappa_0, \\nu_0, \\mathbf{S}_0$ (scale matrix for precision $\\boldsymbol{\\Lambda}$).

Posterior hyperparameters for $(\\boldsymbol{\\mu}_k, \\boldsymbol{\\Lambda}_k = \\boldsymbol{\\Sigma}_k^{-1})$:
$$\\kappa_n = \\kappa_0 + n_k$$
$$\\boldsymbol{m}_n = \\frac{\\kappa_0 \\boldsymbol{m}_0 + n_k \\bar{\\mathbf{y}}_k}{\\kappa_n}$$
$$\\nu_n = \\nu_0 + n_k$$
$$\\mathbf{S}_n = \\mathbf{S}_0 + \\mathbf{S}_{scat,k} + \\frac{\\kappa_0 n_k}{\\kappa_0 + n_k}(\\bar{\\mathbf{y}}_k - \\boldsymbol{m}_0)(\\bar{\\mathbf{y}}_k - \\boldsymbol{m}_0)^T$$

Sampling procedure:
1. Draw precision $\\boldsymbol{\\Lambda}_k | D_k, G_0 \\sim \\text{Wishart}(\\nu_n, \\mathbf{S}_n)$.
2. Draw mean $\\boldsymbol{\\mu}_k | \\boldsymbol{\\Lambda}_k, D_k, G_0 \\sim N_d(\\boldsymbol{m}_n, (\\kappa_n \\boldsymbol{\\Lambda}_k)^{-1})$.")

  sink()
  cat("Generated:", filepath_mvn_post, "\n")
}
cat("DEBUG: Defining function document_predictive_likelihood_math...\n")
# 5. Predictive Likelihood Functions
document_predictive_likelihood_math <- function(base_dir) {
  cat("DEBUG: Executing document_predictive_likelihood_math...\n")
  # Predictive for Gaussian (Normal-Inverse Gamma G0)
  filepath_pred_norm <- file.path(base_dir, "predictive_gaussian_norm_invgamma.md")
  cat("DEBUG: Writing to", filepath_pred_norm, "\n")
  sink(filepath_pred_norm)

  writeLines("# Predictive Likelihood: New Point from Gaussian with Normal-Inverse Gamma Prior

The predictive likelihood for a new data point $y_{new}$, given $G_0 \\sim \\text{Normal-Inverse Gamma}(\\mu_0, \\kappa_0, a_0, b_0)$, is the density of $y_{new}$ under a non-standardized Student's t-distribution:

$$p(y_{new} | G_0) = \\int N(y_{new} | \\mu, \\sigma^2) p(\\mu, \\sigma^2 | G_0) d\\mu d\\sigma^2$$
This evaluates to a density proportional to:
$$p(y_{new} | G_0) \\propto \\left(1 + \\frac{\\kappa_0(y_{new}-\\mu_0)^2}{(\\kappa_0+1)2b_0}\\right)^{-(a_0+1/2)}$$
The exact form used in `Predictive.normal` involves calculating posterior parameters $(a_n, b_n, \\kappa_n)$ as if $y_{new}$ was observed, and then using a ratio of Gamma functions and other terms:
$$p(y_{new}|G_0) = \\frac{\\Gamma(a_n)}{\\Gamma(a_0)} \\frac{b_0^{a_0}}{b_n^{a_n}} \\sqrt{\\frac{\\kappa_0}{\\kappa_n (2\\pi)}}$$ (The $2\\pi$ factor might be handled differently or cancel in practice within the CRP probabilities if not careful with exact density vs. proportional to).
The `dirichletprocess` R code for `Predictive.normal` is: $\\frac{\\Gamma(a_n)}{\\Gamma(a_0)} \\frac{b_0^{a_0}}{b_n^{a_n}} \\sqrt{\\frac{\\kappa_0}{\\kappa_n}}$, where $a_n, b_n, \\kappa_n$ are posterior parameters after observing $y_{new}$. This form is common for the marginal likelihood used in Bayes Factors or model comparison, and is proportional to the Student-t PDF.")

  sink()
  cat("Generated:", filepath_pred_norm, "\n")

  # Predictive for Exponential (Gamma G0)
  filepath_pred_exp <- file.path(base_dir, "predictive_exponential_gamma.md")
  cat("DEBUG: Writing to", filepath_pred_exp, "\n")
  sink(filepath_pred_exp)

  writeLines("# Predictive Likelihood: New Point from Exponential with Gamma Prior

For an Exponential likelihood $F(y | \\lambda_{rate})$ and a Gamma prior $G_0 = \\text{Gamma}(\\lambda_{rate} | a_0, b_0)$ on the rate.

The predictive likelihood for a new data point $y_{new}$ is:
$$p(y_{new} | G_0) = \\int_0^\\infty (\\lambda_{rate} e^{-\\lambda_{rate} y_{new}}) \\frac{b_0^{a_0}}{\\Gamma(a_0)} \\lambda_{rate}^{a_0-1} e^{-b_0 \\lambda_{rate}} d\\lambda_{rate}$$
$$= \\frac{b_0^{a_0}}{\\Gamma(a_0)} \\int_0^\\infty \\lambda_{rate}^{a_0} e^{-(b_0+y_{new})\\lambda_{rate}} d\\lambda_{rate}$$
This is related to the Gamma function. The integral is $\\frac{\\Gamma(a_0+1)}{(b_0+y_{new})^{a_0+1}}$.
So, $$p(y_{new} | G_0) = \\frac{b_0^{a_0}}{\\Gamma(a_0)} \\frac{\\Gamma(a_0+1)}{(b_0+y_{new})^{a_0+1}} = \\frac{a_0 b_0^{a_0}}{(b_0+y_{new})^{a_0+1}}$$
This is the PDF of a Lomax (Pareto Type II) distribution. The R code `Predictive.exponential` calculates $a_n = a_0+1, b_n = b_0+y_{new}$ and then uses $\\frac{\\Gamma(a_n)}{\\Gamma(a_0)} \\frac{b_0^{a_0}}{b_n^{a_n}}$.")

  sink()
  cat("Generated:", filepath_pred_exp, "\n")

  # Predictive for MVN (Normal-Wishart G0)
  filepath_pred_mvn <- file.path(base_dir, "predictive_mvnormal_norm_wishart.md")
  cat("DEBUG: Writing to", filepath_pred_mvn, "\n")
  sink(filepath_pred_mvn)

  writeLines("# Predictive Likelihood: New Point from MVN with Normal-Wishart Prior

For an MVN likelihood $N_d(\\mathbf{y} | \\boldsymbol{\\mu}, \\boldsymbol{\\Sigma})$ and a Normal-Wishart prior $G_0$ on $(\\boldsymbol{\\mu}, \\boldsymbol{\\Lambda}=\\boldsymbol{\\Sigma}^{-1})$ with hyperparameters $(\\boldsymbol{m}_0, \\kappa_0, \\nu_0, \\mathbf{S}_0)$.

The predictive likelihood for a new data point $\\mathbf{y}_{new}$ is the density of $\\mathbf{y}_{new}$ under a multivariate Student's t-distribution:
$$p(\\mathbf{y}_{new} | G_0) = \\int N_d(\\mathbf{y}_{new} | \\boldsymbol{\\mu}, \\boldsymbol{\\Lambda}^{-1}) p(\\boldsymbol{\\mu}, \\boldsymbol{\\Lambda} | G_0) d\\boldsymbol{\\mu} d\\boldsymbol{\\Lambda}$$
This evaluates to $t_{\\nu_0-d+1}\\left(\\mathbf{y}_{new} \\Big| \\boldsymbol{m}_0, \\frac{\\mathbf{S}_0^{-1}(\\kappa_0+1)}{\\kappa_0(\\nu_0-d+1)}\\right)$ (scaled and shifted multivariate t).

The R code `Predictive.mvnormal` calculates posterior parameters $(\\boldsymbol{m}_n, \\kappa_n, \\nu_n, \\mathbf{S}_n)$ as if $\\mathbf{y}_{new}$ were observed, and then uses a formula involving ratios of determinants and Gamma functions:
$$p(\\mathbf{y}_{new} | G_0) = \\pi^{-d/2} \\frac{\\Gamma_d(\\nu_n/2)}{\\Gamma_d(\\nu_0/2)} \\frac{|\\mathbf{S}_0|^{\\nu_0/2}}{|\\mathbf{S}_n|^{\\nu_n/2}} \\left(\\frac{\\kappa_0}{\\kappa_n}\\right)^{d/2}$$
where $\\Gamma_d(\\cdot)$ is the multivariate Gamma function. This is proportional to the multivariate Student-t PDF.")

  sink()
  cat("Generated:", filepath_pred_mvn, "\n")
}
cat("DEBUG: Defining function document_mcmc_components_math...\n")
# 6. MCMC Components (CRP, MH)
document_mcmc_components_math <- function(base_dir) {
  cat("DEBUG: Executing document_mcmc_components_math...\n")
  # CRP details (already in mathematical_basis.md, but can be expanded here if needed)
  # Metropolis-Hastings
  filepath_mh <- file.path(base_dir, "metropolis_hastings.md")
  cat("DEBUG: Writing to", filepath_mh, "\n")
  sink(filepath_mh)

  writeLines("# Metropolis-Hastings Algorithm for Non-Conjugate Parameter Updates

When updating cluster parameters $\\phi_k$ for a non-conjugate model, we sample from the posterior $p(\\phi_k | D_k, G_0) \\propto \\left( \\prod_{y_j \\in D_k} F(y_j | \\phi_k) \\right) p_0(\\phi_k)$, where $p_0(\\phi_k)$ is the prior density from $G_0$.

Algorithm steps (for one cluster $k$ at MCMC iteration $t+1$):
1. Current state: $\\phi_k^{(t)}$.
2. Propose a new state $\\phi_k'$ from a proposal distribution $q(\\phi_k' | \\phi_k^{(t)})$. A common choice is a symmetric random walk, e.g., $\\phi_k' \\sim N(\\phi_k^{(t)}, \\Sigma_{prop})$.
3. Calculate the acceptance ratio $\\mathcal{A}(\\phi_k', \\phi_k^{(t)})$:
   $$\\mathcal{A}(\\phi_k', \\phi_k^{(t)}) = \\min\\left(1, \\frac{p(\\phi_k' | D_k, G_0) q(\\phi_k^{(t)} | \\phi_k')}{p(\\phi_k^{(t)} | D_k, G_0) q(\\phi_k' | \\phi_k^{(t)})}\\right)$$
   If $q$ is symmetric (e.g., Gaussian random walk), $q(\\phi_k^{(t)} | \\phi_k') = q(\\phi_k' | \\phi_k^{(t)})$, so these terms cancel:
   $$\\mathcal{A}(\\phi_k', \\phi_k^{(t)}) = \\min\\left(1, \\frac{\\left( \\prod_{y_j \\in D_k} F(y_j | \\phi_k') \\right) p_0(\\phi_k')}{\\left( \\prod_{y_j \\in D_k} F(y_j | \\phi_k^{(t)}) \\right) p_0(\\phi_k^{(t)})}\\right)$$
4. Draw $u \\sim U(0,1)$.
5. If $u < \\mathcal{A}(\\phi_k', \\phi_k^{(t)})$, set $\\phi_k^{(t+1)} = \\phi_k'$ (accept).
6. Else, set $\\phi_k^{(t+1)} = \\phi_k^{(t)}$ (reject).

Implementation typically uses log-probabilities for numerical stability.")

  sink()
  cat("Generated:", filepath_mh, "\n")
}
cat("DEBUG: Defining function document_concentration_param_math...\n")
# 7. Concentration Parameter Update
document_concentration_param_math <- function(filepath) {
  cat("DEBUG: Executing document_concentration_param_math for", filepath, "\n")
  sink(filepath)

  writeLines("# Concentration Parameter $\\alpha$ Update (Escobar & West, 1995)

The concentration parameter $\\alpha$ controls the expected number of clusters. It is typically assigned a Gamma prior: $\\alpha \\sim \\text{Gamma}(a_\\alpha, b_\\alpha)$.

Given $K$ current clusters among $N$ data points, the update uses an auxiliary variable $\\eta$. Following Escobar and West (1995), and as implemented in `dirichletprocess` (see `update_concentration` R function):

1. Draw $\\eta | \\alpha, N \\sim \\text{Beta}(\\alpha + 1, N)$.
2. The posterior for $\\alpha$ is then a mixture of two Gamma distributions:
   $$p(\\alpha | K, N, \\eta, a_\\alpha, b_\\alpha) = \\pi_\\eta \\cdot \\text{Gamma}(\\alpha | a_\\alpha+K, b_\\alpha-\\log(\\eta)) + (1-\\pi_\\eta) \\cdot \\text{Gamma}(\\alpha | a_\\alpha+K-1, b_\\alpha-\\log(\\eta))$$
   where the mixture probability $\\pi_\\eta$ is given by:
   $$\\pi_\\eta = \\frac{a_\\alpha+K-1}{(a_\\alpha+K-1) + N(b_\\alpha-\\log(\\eta))}$$
   (This assumes $a_\\alpha+K-1 > 0$. If $K=0$, the second Gamma component might use shape $a_\\alpha$. The R code handles cases where $a_\\alpha+K-1 \\le 0$ by adjusting which component is chosen or falling back to the prior).

A draw is made from this mixture distribution to get the new $\\alpha$.")

  sink()
  cat("Generated:", filepath, "\n")
}
cat("DEBUG: Defining function document_hierarchical_math...\n")
# 8. Hierarchical Dirichlet Process Specific Math
document_hierarchical_math <- function(base_dir) {
  cat("DEBUG: Executing document_hierarchical_math...\n")
  filepath_hdp <- file.path(base_dir, "hierarchical_dp_math.md")
  cat("DEBUG: Writing to", filepath_hdp, "\n")
  sink(filepath_hdp)

  writeLines("# Hierarchical Dirichlet Process (HDP) - Mathematical Details

The HDP (Teh et al., 2005) allows sharing of statistical strength across multiple related groups, where each group $j$ has its own DP mixture model, but these DPs share a common, global base measure $G_0$, which is itself drawn from a DP.

## Model Specification

For group $j=1, \\dots, J$:
$$y_{ji} | \\theta_{ji} \\sim F(\\cdot | \\theta_{ji}) \\quad \\text{(Observation } i \\text{ in group } j)$$
$$\\theta_{ji} | G_j \\sim G_j \\quad \\text{(Parameters for observation } ji)$$
$$G_j | G_0, \\alpha_0 \\sim DP(G_0, \\alpha_0) \\quad \\text{(Group-specific DP, } \\alpha_0 \\text{ is concentration for } G_j)$$
$$G_0 | H, \\gamma \\sim DP(H, \\gamma) \\quad \\text{(Global DP, } \\gamma \\text{ is top-level concentration})$$

Where:
- $F$ is the likelihood function.
- $G_j$ is the random measure for group $j$.
- $G_0$ is the global random measure, serving as the base measure for each $G_j$. It is discrete: $G_0 = \\sum_{k=1}^{\\infty} \\beta_k \\delta_{\\phi_k}$.
- $\\phi_k \\sim H$ are the global atoms (parameters) drawn from a hyperprior $H$.
- $\\boldsymbol{\\beta} = (\\beta_1, \\beta_2, \\dots)$ are the global stick-breaking weights from $DP(H, \\gamma)$, i.e., $\\beta_k = v_k \\prod_{l<k}(1-v_l)$ with $v_k \\sim \\text{Beta}(1, \\gamma)$.
- $\\alpha_0$ is the concentration parameter for group-level DPs (can be different per group, $\\alpha_j$).

## Key Update Steps in HDP Gibbs Sampler (Simplified, based on Neal's methods where applicable)

### 1. Update Data Point Assignments $c_{ji}$ (within each group $j$)
For each data point $y_{ji}$ in group $j$, assign it to a local cluster $l$ (with params $\\theta_{jl}$) or a new local cluster using the CRP specific to $G_j$. The probability of assigning to an existing local cluster $l$ (with $n_{jl,-i}$ points) is proportional to $n_{jl,-i} F(y_{ji} | \\theta_{jl})$. The probability of assigning to a new local cluster is proportional to $\\alpha_0 \\int F(y_{ji} | \\theta) dG_0(\\theta)$. The parameters $\\theta_{jl}$ are themselves draws from $G_0$.

### 2. Update Local Cluster Parameters $\\theta_{jl}$ (within each group $j$)
For each local cluster $l$ in group $j$, its parameters $\\theta_{jl}$ are effectively a draw from the global $G_0$. In practice, these are often linked to one of the global atoms $\\phi_k$. The update involves re-sampling which global atom $\\phi_k$ this local cluster type $\\theta_{jl}$ corresponds to, with probability proportional to $\\beta_k \\prod_{y \\in D_{jl}} F(y | \\phi_k)$.

### 3. Update Global Atoms $\\phi_k$
Each global atom $\\phi_k$ is updated based on all data points $D_k = \\{ y_{ji} : \\theta_{ji} \\text{ is an instance of } \\phi_k \\}$ that are associated with it, across all groups. The update is from the posterior:
$$\\phi_k | D_k, H \\sim p(\\phi_k | D_k, H) \\propto \\left( \\prod_{y \\in D_k} F(y | \\phi_k) \\right) p_H(\\phi_k)$$
where $p_H(\\phi_k)$ is the prior density from $H$. This step uses standard conjugate updates or MH if $H$ is not conjugate to $F$.

### 4. Update Global Stick-Breaking Weights $\\boldsymbol{\\beta}$ (and potentially add/remove atoms from $G_0$)
Let $m_k$ be the number of distinct local cluster types (tables) across all groups that are instances of the global atom $\\phi_k$. Let $K_0$ be the number of currently instantiated global atoms.
The posterior for the weights $(\\beta_1, \\dots, \\beta_{K_0}, \\beta_{new_mass})$ is Dirichlet:
$$(\\beta_1, \\dots, \\beta_{K_0}, \\sum_{l=K_0+1}^{\\infty} \\beta_l) \\sim \\text{Dirichlet}(m_1, \\dots, m_{K_0}, \\gamma)$$
The sum $\\sum_{l=K_0+1}^{\\infty} \\beta_l$ is the mass for new global atoms. New atoms $\\phi_{K_0+1}, \\dots$ can be drawn from $H$, and their weights derived by further stick-breaking on this remaining mass, often with truncation.

### 5. Update Concentration Parameters $\\alpha_0$ (group-level) and $\\gamma$ (global-level)
- Each $\\alpha_j$ (or common $\\alpha_0$) is updated using the Escobar & West method, based on the number of local clusters $K_j$ and data points $N_j$ within group $j$.
- The global $\\gamma$ is updated using Escobar & West, based on the number of distinct global atoms $K_0$ currently used by any table, and the total number of tables $M = \\sum_k m_k$ across all groups.

This is a simplified overview. Specific algorithms (e.g., Neal's Algorithm 8 for HDPs, or direct assignment methods) have more detailed steps for managing assignments to global atoms and sampling new atoms.")

  sink()
  cat("Generated:", filepath_hdp, "\n")
}
# --- Call the main documentation function ---
# To run, uncomment the line below and ensure your working directory allows writing to 'analysis/code_documentation/'
document_all_maths()
cat("R script 'mathematical_documentation_generator.R' is defined.\nCall 'document_all_maths()' to generate the Markdown files.\n")
