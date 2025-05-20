# dirichletprocess Package Analysis

## Package Structure

Package Structure Overview
 =======================
 Total R files:  65 
 Total functions:  161 


## Core Workflow

# dirichletprocess Core Workflow

## 1. Object Creation
The main entry point is typically through one of these constructors:
- `DirichletProcessGaussian()`
- `DirichletProcessBeta()`
- `DirichletProcessMvnormal()`
- `DirichletProcessWeibull()`

These functions all call `DirichletProcessCreate()` which sets up the basic object structure.

## 2. Initialization
The `Initialise()` function prepares the object for fitting:
- For conjugate mixtures: `Initialise.conjugate()`
- For non-conjugate mixtures: `Initialise.nonconjugate()`

## 3. Fitting
The `Fit()` function runs the MCMC sampling:
- For standard DPs: `Fit.default()`
- For hierarchical DPs: `Fit.hierarchical()`

## 4. Core MCMC Components
Each iteration of the MCMC process involves:
- `ClusterComponentUpdate()`: Updates cluster assignments
- `ClusterParameterUpdate()`: Updates cluster parameters
- `UpdateAlpha()`: Updates concentration parameter



## S3 Class Structure

# S3 Class Structure in dirichletprocess

## Main S3 Classes

### dirichletprocess
Base class for all Dirichlet process objects

#### Subclasses based on mixture distribution:
- `normal`: Gaussian mixture
- `beta`: Beta mixture
- `mvnormal`: Multivariate normal mixture
- `mvnormal2`: Alternative multivariate normal implementation
- `weibull`: Weibull mixture
- `exponential`: Exponential mixture

#### Subclasses based on implementation type:
- `conjugate`: For conjugate prior-likelihood pairs
- `nonconjugate`: For non-conjugate prior-likelihood pairs
- `hierarchical`: For hierarchical Dirichlet processes
- `markov`: For hidden Markov models

### MixingDistribution
Class for mixture distribution objects

## Method Dispatch

The package uses S3 method dispatch extensively. Key generics include:

### Fit
Methods:
- `Fit.default`
- `Fit.hierarchical`
- `Fit.markov`

### ClusterComponentUpdate
Methods:
- `ClusterComponentUpdate.conjugate`
- `ClusterComponentUpdate.hierarchical`
- `ClusterComponentUpdate.nonconjugate`

### ClusterParameterUpdate
Methods:
- `ClusterParameterUpdate.conjugate`
- `ClusterParameterUpdate.nonconjugate`

### PriorDraw
Methods:
- `PriorDraw.beta`
- `PriorDraw.beta2`
- `PriorDraw.exponential`
- `PriorDraw.hierarchical`
- `PriorDraw.mvnormal`
- `PriorDraw.mvnormal2`
- `PriorDraw.normal`
- `PriorDraw.normalFixedVariance`
- `PriorDraw.weibull`

### PosteriorDraw
Methods:
- `PosteriorDraw.exponential`
- `PosteriorDraw.mvnormal`
- `PosteriorDraw.mvnormal2`
- `PosteriorDraw.nonconjugate`
- `PosteriorDraw.normal`
- `PosteriorDraw.normalFixedVariance`
- `PosteriorDraw.weibull`

### UpdateAlpha
Methods:
- `UpdateAlpha.default`
- `UpdateAlpha.hierarchical`

### Likelihood
Methods:
- `Likelihood.beta`
- `Likelihood.beta2`
- `Likelihood.exponential`
- `Likelihood.mvnormal`
- `Likelihood.mvnormal2`
- `Likelihood.normal`
- `Likelihood.normalFixedVariance`
- `Likelihood.weibull`



## Mathematical Foundation

# Mathematical Foundation of dirichletprocess Algorithms

## Dirichlet Process Mixture Model (DPMM)

The model assumes data $y_i$ are generated from a mixture model where the mixing distribution $G$ itself is drawn from a Dirichlet Process (DP). The hierarchical representation is:

$$y_i | \theta_i \sim F(\cdot | \theta_i)$$
$$\theta_i | G \sim G$$
$$G | G_0, \alpha \sim DP(G_0, \alpha)$$

Where:
- $y_i$ is the $i$-th observation.
- $F(\cdot | \theta_i)$ is the likelihood function for an observation given parameters $\theta_i$ (e.g., Gaussian $N(y_i | \mu_i, \sigma_i^2)$, Beta $Beta(y_i | a_i, b_i)$).
- $\theta_i$ are the parameters specific to the $i$-th observation.
- $G$ is a random probability measure (the mixing distribution) drawn from a Dirichlet Process.
- $G_0$ is the base measure, representing the prior belief about the distribution of parameters $\theta$. It is the expected value of $G$, $E[G] = G_0$.
- $\alpha > 0$ is the concentration parameter. Larger $\alpha$ leads to $G$ being closer to $G_0$ and implies a larger number of distinct parameter values (clusters) a priori.

Due to the discreteness of draws from a DP, many $\theta_i$ will share common values. Let $\{\phi_k\}_{k=1}^K$ be the set of unique parameter values, where $K$ is the number of clusters.

## Chinese Restaurant Process (CRP) / Gibbs Sampling for Cluster Assignments (Neal's Algorithm 3)

The Chinese Restaurant Process (CRP) formulation provides a constructive way to understand the clustering property of the DP. For Gibbs sampling, we integrate out $G$. The conditional posterior probability for assigning data point $y_i$ to cluster $k$ (parameterized by $\phi_k$), given all other assignments $c_{-i}$ and current cluster parameters, is (Neal, 2000, Algorithm 3):

$$P(c_i = k | c_{-i}, \mathbf{y}, \alpha, G_0, \{\phi_j\}_{j \neq k}) \propto \begin{cases}
n_{k,-i} \cdot F(y_i | \phi_k) & \text{for an existing cluster } k \text{ (where } n_{k,-i} > 0) \\
\alpha \cdot \int F(y_i | \phi) dG_0(\phi) & \text{for a new cluster (i.e., } k = K_{new})
\end{cases}$$

Where:
- $c_i$ is the cluster assignment for data point $y_i$.
- $c_{-i}$ denotes the cluster assignments for all data points except $y_i$.
- $n_{k,-i}$ is the number of data points in cluster $k$, excluding point $y_i$.
- $F(y_i | \phi_k)$ is the likelihood of data point $y_i$ given the parameters $\phi_k$ of cluster $k$.
- $\alpha$ is the concentration parameter.
- $\int F(y_i | \phi) dG_0(\phi)$ is the marginal likelihood of $y_i$ under the base measure $G_0$, often called the predictive likelihood. This is used when proposing a new cluster, whose parameters $\phi_{K_{new}}$ would be drawn from $G_0$ (or a posterior based on $y_i$ if conjugate).

This is Neal's (2000) Algorithm 3. Algorithm 2 is similar but marginalizes over $\phi_k$ as well, which is common for fully conjugate models.

## Cluster Parameter Updates (Neal's Algorithm 3)

After re-assigning all $c_i$, the parameters $\phi_k$ for each cluster $k$ (now containing points $D_k = \{y_i : c_i = k\}$) are updated.

### Conjugate Case

If $G_0$ is conjugate to the likelihood $F$, the parameters $\phi_k$ for each cluster $k$ are drawn from their posterior distribution:

$$\phi_k | D_k, G_0 \sim p(\phi_k | D_k, G_0)$$

This posterior is analytically known. For example, if $F$ is Gaussian and $G_0$ is Normal-Inverse-Gamma, the posterior $p(\phi_k | D_k, G_0)$ is also Normal-Inverse-Gamma with updated hyperparameters.

### Non-conjugate Case (Neal's Algorithm 8 for assignments, Metropolis-Hastings for parameters)

When $G_0$ is not conjugate to $F$, sampling $\phi_k$ directly from the posterior is not possible. Instead, MCMC methods like Metropolis-Hastings are used. The target density is:

$$p(\phi_k | D_k, G_0) \propto \left( \prod_{y_j \in D_k} F(y_j | \phi_k) \right) p_0(\phi_k)$$

where $p_0(\phi_k)$ is the prior density of $\phi_k$ from $G_0$.
For cluster assignments in non-conjugate cases, Neal's Algorithm 8 is often used, which introduces $m$ auxiliary parameters drawn from $G_0$ to help propose new clusters, avoiding the direct calculation of $\int F(y_i | \phi) dG_0(\phi)$.

## Concentration Parameter $\alpha$ Update (Escobar & West, 1995)

The concentration parameter $\alpha$ can be given a prior (e.g., Gamma distribution $p(\alpha) = \text{Gamma}(\alpha | a_\alpha, b_\alpha)$) and updated using an auxiliary variable method (Escobar and West, 1995). Given $K$ current clusters and $N$ data points, the update involves drawing an auxiliary variable $\eta \sim \text{Beta}(\alpha+1, N)$. The posterior for $\alpha$ is then a mixture of two Gamma distributions:

$$p(\alpha | K, N, \eta, \text{priors}) = \pi_\eta \text{Gamma}(\alpha | a_\alpha+K, b_\alpha-\log(\eta)) + (1-\pi_\eta) \text{Gamma}(\alpha | a_\alpha+K-1, b_\alpha-\log(\eta))$$

Where the mixture weight $\pi_\eta$ depends on $a_\alpha, K, N, b_\alpha,$ and $\log(\eta)$.



## Performance Analysis

# dirichletprocess Package: Profiling Results and Bottleneck Analysis

**Author**: Analysis based on profiling scripts
**Date**: May 19, 2025

## 1. Introduction

This report summarizes the performance profiling results for the `dirichletprocess` R package. The goal of this profiling was to identify computational bottlenecks within the core MCMC sampling algorithms, specifically to guide the C++ reimplementation efforts. Profiling was conducted using `profvis` for overall function tracing and `microbenchmark` for timing specific components. Different mixture types (Gaussian, Beta, MVN) and varying data sizes were considered.

The key outputs analyzed are:

* `benchmark_summary.txt`: Contains `microbenchmark` timings for `ClusterComponentUpdate`, `ClusterParameterUpdate`, and `UpdateAlpha`.

* `profvis` HTML outputs (e.g., `gaussian_n100.html`, `beta_n100.html`, etc.): Provide detailed flame graphs and data views of the `Fit` function execution.

## 2. Microbenchmark Component Analysis

The `microbenchmark` results provide precise timings for the core iterative components of the MCMC sampler.

**(Note: The following section reproduces the key insights from the provided `benchmark_summary.txt`)**

```text
# Benchmark Results


##  Gaussian ClusterComponentUpdate (n=200) 
                                       expr    min     lq     mean
1 ClusterComponentUpdate(dp_gaussian_bench) 9.0073 9.1943 9.592545
   median     uq     max neval
1 9.33175 9.4064 14.6837    20


##  Gaussian ClusterParameterUpdate (n=200) 
                                       expr   min     lq    mean
1 ClusterParameterUpdate(dp_gaussian_bench) 106.1 107.75 117.405
  median    uq   max neval
1 108.65 113.1 255.4    20


##  Gaussian UpdateAlpha (n=200) 
                            expr  min   lq  mean median   uq   max
1 UpdateAlpha(dp_gaussian_bench) 16.5 16.7 19.76     17 17.4 106.7
  neval
1    50


##  Beta ClusterComponentUpdate (n=100, Non-conjugate) 
                                   expr    min     lq    mean
1 ClusterComponentUpdate(dp_beta_bench) 7.2424 7.2891 7.41192
  median     uq    max neval
1 7.3018 7.4016 8.1218    10


##  Beta ClusterParameterUpdate (n=100, Non-conjugate) 
                                   expr    min     lq     mean
1 ClusterParameterUpdate(dp_beta_bench) 488.15 489.75 498.9444
  median     uq    max neval
1 491.40 509.00 516.21     5


##  Beta UpdateAlpha (n=100, Non-conjugate) 
                            expr    min     lq   mean median    uq
1 UpdateAlpha(dp_beta_bench) 18.1 18.35 20.494  18.85 19.55
    max neval
1 113.4    50
```

### Key Findings from Microbenchmarks:

#### Gaussian Mixture (Conjugate, n=200)

* **`UpdateAlpha`**: Median ~17.0 µs. Very fast, as expected.

* **`ClusterComponentUpdate`**: Median ~9.33 ms. This is a significant portion of the work, involving iteration over data points and clusters.

* **`ClusterParameterUpdate`**: Median ~108.65 ms. This is the most time-consuming component for the conjugate Gaussian model, involving posterior draws for each active cluster.

#### Beta Mixture (Non-conjugate, n=100)

* **`UpdateAlpha`**: Median ~18.85 µs. Remains very fast.

* **`ClusterComponentUpdate`**: Median ~7.30 ms. Comparable in impact to the Gaussian case, considering the smaller `n`.

* **`ClusterParameterUpdate`**: Median ~491.40 ms. This is **extremely slow**, nearly 0.5 seconds per call for only 100 data points. It's roughly 4.5 times slower than its conjugate Gaussian counterpart, even with half the data size.

**Conclusion from Microbenchmarks:**
The `ClusterParameterUpdate` step is the most computationally intensive, especially for non-conjugate mixtures like Beta, where it's an order of magnitude slower. `ClusterComponentUpdate` is the second major contributor to runtime. `UpdateAlpha` is negligible in comparison.

## 3. `profvis` Full Sampler Analysis (Inferred)

The `profvis` HTML outputs (e.g., `gaussian_n1000.html`, `beta_n500.html`) provide a holistic view of the `Fit` function. While not directly embedded here, common patterns emerge:

* The vast majority of the execution time within `Fit()` is spent inside its main MCMC loop.

* Within this loop, the calls to `ClusterParameterUpdate()` and `ClusterComponentUpdate()` (and their S3 dispatched methods) will dominate the flame graph.

### Conjugate Models (Gaussian, MVN)

* **`ClusterParameterUpdate.*`**:

  * `PosteriorDraw.*`: Calculating posterior hyperparameters (which involves iterating/summing over data points within each cluster) and then sampling from the known posterior distribution. For MVN models, matrix operations within `PosteriorDraw.mvnormal` (e.g., sums of squares and cross-products, inversions if not careful, Cholesky decompositions) would be significant.

* **`ClusterComponentUpdate.*`**:

  * `Likelihood.*`: These functions are called $N \times K_{avg}$ times per MCMC iteration (where $N$ is the number of data points, $K_{avg}$ is the average number of clusters). For MVN, `mvtnorm::dmvnorm` is a key function.

  * `Predictive.*`: Called $N$ times per MCMC iteration for evaluating new cluster probabilities.

  * `ClusterLabelChange.*`: Management of R lists for `clusterParameters` and `pointsPerCluster` (especially resizing and subsetting) can introduce overhead, particularly if the number of clusters changes frequently.

### Non-Conjugate Models (Beta)

* **`ClusterParameterUpdate.nonconjugate`**: This is the **overwhelming bottleneck**.

  * The primary reason is the call to `PosteriorDraw.nonconjugate`, which internally uses `MetropolisHastings`.

  * **`MetropolisHastings` Loop**: This loop runs for `dpObj$mhDraws` (e.g., 250) iterations *for each active cluster* in *each MCMC iteration of `Fit`*.

    * **`Likelihood.beta`**: Critically, this is called for *all data points assigned to the current cluster* at *every single Metropolis-Hastings step*. Summing log-likelihoods across these points is the most intensive part.

    * `PriorDensity.beta`: Called once per MH step.

    * `MhParameterProposal.beta`: Called once per MH step.

  * **Low Metropolis-Hastings Acceptance Ratios**: The console output showed acceptance ratios around 0.3% - 0.9% for Beta mixtures. This means over 99% of the computationally expensive MH proposals (including likelihood calculations) are rejected. This highlights an inefficiency in the MH sampler's proposal/tuning for the Beta model, making the already expensive parameter update step even slower in terms of effective samples.

* **`ClusterComponentUpdate.nonconjugate`**:

  * `Likelihood.beta`: Still called frequently ($N \times K_{avg}$ for existing clusters, and $N \times m_{aux}$ for auxiliary parameters).

  * `PriorDraw.beta`: Used for drawing auxiliary parameters.

  * `ClusterLabelChange.nonconjugate`.

## 4. Summary of Identified Bottlenecks for C++ Optimization

Based on the microbenchmarks and the expected behavior from `profvis` analysis, the following are the primary targets for C++ optimization, in approximate order of impact:

1. **`ClusterParameterUpdate` (Non-Conjugate Models)**:

   * The entire **`MetropolisHastings` sampling loop** within `PosteriorDraw.nonconjugate`.

   * Specifically, the repeated calculation of **`Likelihood.*`** (e.g., `Likelihood.beta`) over all data points in a cluster, for each MH step.

   * Also, `PriorDensity.*` and `MhParameterProposal.*` within the MH loop.

   * *Goal*: Make each MH step significantly faster.

2. **`Likelihood.*` Functions (All Kernels)**:

   * These are the most frequently called low-level computational routines.

   * Optimizing these in C++ will benefit both `ClusterComponentUpdate` and `ClusterParameterUpdate` (especially non-conjugate).

   * Focus on efficient mathematical computations and vectorized operations where possible.

3. **`ClusterComponentUpdate` Core Logic (All Models)**:

   * The outer loop over $N$ data points.

   * The inner loop over $K$ clusters for calculating probabilities (heavy use of `Likelihood`).

   * Calculation of new cluster probabilities (using `Predictive` for conjugate, or `Likelihood` with auxiliary parameters for non-conjugate).

   * The categorical sampling step.

4. **`PosteriorDraw.*` and `Predictive.*` (Conjugate Models)**:

   * The calculation of posterior hyperparameters (often involves sums/means over data in clusters).

   * The actual random sampling from these posteriors.

   * Mathematical formulas for predictive distributions.

   * For MVN, efficient matrix algebra is key.

5. **`ClusterLabelChange.*`**:

   * Efficient management of R data structures (cluster parameters, labels, counts) when clusters are added, removed, or re-indexed. C++ data structures like `std::vector` can offer better performance for dynamic resizing and access compared to R lists if large parts of this logic are moved to C++.

6. **Data Structure Access and Management**:

   * Minimizing copying of data between R and C++.

   * Using efficient C++ data structures internally.

   * Reducing overhead from repeated S3 method dispatch for core computational routines by having C++ handle larger chunks of work.

## 5. Conclusion

The profiling results clearly indicate that the `ClusterParameterUpdate` step, particularly for non-conjugate mixtures due to the nested Metropolis-Hastings sampling, is the most significant performance bottleneck. The `Likelihood` functions, being called extensively within both `ClusterParameterUpdate` and `ClusterComponentUpdate`, are also critical targets.

Reimplementing these core components and their dependencies (like likelihood, prior, and posterior draw functions) in C++ is expected to yield substantial performance improvements. Additionally, addressing the low MH acceptance ratio for the Beta mixture through better proposal tuning (an R-level concern) would complement the C++ speedups for practical usability.


## Optimization Plan

# C++ Optimization Roadmap

## Priority 1: Core Data Structures

1. DP Object representation
2. Mixing Distribution representation
3. Cluster storage structures

## Priority 2: Likelihood Functions

1. Gaussian likelihood
2. Multivariate Normal likelihood
3. Beta likelihood
4. Weibull likelihood
5. Exponential likelihood

## Priority 3: Core MCMC Algorithms

1. ClusterComponentUpdate.conjugate
2. ClusterComponentUpdate.nonconjugate
3. ClusterParameterUpdate.conjugate
4. ClusterParameterUpdate.nonconjugate
5. UpdateAlpha

## Priority 4: Prior and Posterior Sampling

1. PriorDraw functions
2. PosteriorDraw functions
3. Metropolis-Hastings implementation

## Priority 5: Integration

1. R/C++ interface functions
2. S3 dispatch bridging
3. Error handling and validation

## Expected Performance Gains

1. Likelihood calculations: 50-100x speedup
2. ClusterComponentUpdate: 20-50x speedup
3. ClusterParameterUpdate: 10-30x speedup
4. Overall MCMC: 15-40x speedup

## Implementation Timeline

1. Weeks 1-2: Core data structures and likelihood functions
2. Weeks 3-4: ClusterComponentUpdate for conjugate models
3. Weeks 5-6: ClusterParameterUpdate and non-conjugate implementations
4. Weeks 7-8: Integration and testing
