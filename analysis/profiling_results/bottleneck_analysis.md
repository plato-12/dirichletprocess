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
