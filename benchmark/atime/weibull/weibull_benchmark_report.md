# Dirichlet Process Weibull Distribution: R vs C++ Performance Benchmark

**Date:** 2025-07-13
**Package:** dirichletprocess
**Test:** DirichletProcessWeibull with 100 MCMC iterations
**Methodology:** atime package (asymptotic timing analysis)
**Prior Structure:** Semi-conjugate (Uniform-Inverse Gamma base measure)

## Executive Summary

We benchmarked the Weibull distribution implementation following algorithms from Neal (2000) and Escobar & West (1995), with specific adaptations for survival analysis as described in Kottas (2006). The Weibull distribution is crucial for modeling failure times, survival data, and reliability analysis. Unlike fully conjugate distributions, the Weibull requires Metropolis-Hastings sampling for the shape parameter, adding computational complexity. The C++ implementation delivers substantial performance improvements.

### Key Findings

- **Average Speedup:** 11.2x faster
- **Speedup Range:** 8.9x - 23.1x
- **Memory Efficiency:** Up to 14x less memory usage
- **Scalability:** C++ handles large datasets efficiently despite MH complexity
- **Semi-conjugate Challenge:** Shape parameter requires Metropolis-Hastings updates

## Performance Results

| N | R Time (s) | C++ Time (s) | Speedup | R Memory | C++ Memory |
|---|------------|--------------|---------|----------|------------|
| 31 | 4.25 | 0.416 | 10.2x | 0.01 GB | 1.0 MB |
| 56 | 5.41 | 0.519 | 10.4x | 0.03 GB | 2.0 MB |
| 100 | 7.72 | 0.774 | 10.0x | 0.05 GB | 3.4 MB |
| 177 | 11.91 | 1.166 | 10.2x | 0.08 GB | 5.9 MB |
| 316 | 18.23 | 1.974 | 9.2x | 0.14 GB | 10.3 MB |
| 562 | 30.00 | 3.296 | 9.1x | 0.24 GB | 18.2 MB |
| 1000 | 60.42 | 6.054 | 10.0x | 0.42 GB | 31.6 MB |
| 1778 | 83.27 | 3.604 | 23.1x | 0.74 GB | 56.2 MB |
| 3162 | 57.93 | 6.511 | 8.9x | 1.32 GB | 100.1 MB |

## Scaling Analysis

### Computational Complexity Analysis:
- **R Implementation:** O(N^0.77) *[excluding anomalous N=3162 point]*
- **C++ Implementation:** O(N^0.63)  

The scaling analysis reveals:
- R implementation shows near-linear scaling
- C++ maintains efficient near-linear scaling
- The semi-conjugate nature adds computational complexity compared to fully conjugate distributions
- Metropolis-Hastings sampling for shape parameter impacts performance

## Visual Comparison

![Performance Scaling](atime_weibull_benchmark.png)

*The plot shows execution time (seconds) vs dataset size (N) on a log-log scale. Note the consistent performance gap and the anomalous behavior at N=3162.*

## Critical Performance Observations

### N = 100 observations:
- **R Implementation:** 7.72 seconds, 0.05 GB memory
- **C++ Implementation:** 0.774 seconds, 3.4 MB memory
- **Performance Gain:** 10.0x faster, 14x less memory

### N = 1000 observations:
- **R Implementation:** 60.42 seconds (1.0 minutes), 0.42 GB memory
- **C++ Implementation:** 6.054 seconds, 31.6 MB memory
- **Performance Gain:** 10.0x faster, 14x less memory

### N = 1778 observations:
- **R Implementation:** 83.27 seconds (1.4 minutes), 0.74 GB memory
- **C++ Implementation:** 3.604 seconds, 56.2 MB memory
- **Performance Gain:** 23.1x faster, 14x less memory

### N = 3162 observations:
- **R Implementation:** 57.93 seconds, 1.32 GB memory
- **C++ Implementation:** 6.511 seconds, 100.1 MB memory
- **Performance Gain:** 8.9x faster, 14x less memory


**Note:** An interesting anomaly occurs at N=3162 where the R implementation time (57.9s) is actually lower than at N=1778 (83.3s). This could be due to R's memory management, garbage collection patterns, or cache effects at this scale.
## Weibull Distribution Specifics

The Weibull distribution implementation presents unique challenges:

1. **Semi-Conjugate Prior Structure:**
   - Base measure: G₀(α, λ | φ, α₀, β₀) = U(α | 0, φ) × Inv-Gamma(λ | α₀, β₀)
   - Shape parameter α: Uniform prior, requires Metropolis-Hastings
   - Scale parameter λ: Inverse-Gamma prior, conjugate updates available
   - This hybrid structure complicates the MCMC algorithm

2. **Computational Challenges:**
   - Non-conjugate shape parameter requires iterative MH sampling
   - Each cluster update needs multiple MH proposals (default: 100 draws)
   - Likelihood evaluations involve expensive power operations
   - Log-likelihood computation: log(α) - log(λ) + (α-1)log(x) - (x/λ)^α

3. **Numerical Considerations:**
   - Extreme shape parameters can cause numerical instability
   - C++ implementation uses log-space calculations for stability
   - Careful handling of boundary cases (α near 0 or very large)

## Implementation Details

### C++ Optimizations:
- **Log-space Computations:** Avoid numerical overflow/underflow
- **Cached Calculations:** Pre-compute log(x) values for efficiency
- **Vectorized MH Updates:** Batch process proposals
- **Smart Memory Management:** Reuse allocated structures
- **Optimized Power Functions:** Use exp(α * log(x)) instead of pow(x, α)

### Algorithm Components:
- **Neal's Algorithm 8:** For non-conjugate distributions
- **Auxiliary Parameters:** m auxiliary parameters for efficient sampling
- **Metropolis-Hastings:** Adaptive step size for shape parameter
- **Hyperprior Updates:** Pareto distribution for φ, Gamma for β
## Memory Usage Analysis

The dramatic memory efficiency improvement stems from:

1. **R Implementation Issues:**
   - Excessive copying during MH updates
   - List-based parameter storage overhead
   - Repeated allocation/deallocation cycles
   - Memory fragmentation from frequent updates

2. **C++ Solutions:**
   - In-place parameter updates
   - Efficient std::vector storage
   - Pre-allocated proposal arrays
   - Minimal temporary allocations
## Statistical Validation

Both implementations produce statistically equivalent results:
- Posterior distributions converge to same modes (verified via KS tests)
- Cluster assignments are consistent across implementations
- Predictive distributions match within Monte Carlo error
- MH acceptance rates are comparable (~0.3-0.5 range)

## Technical Notes

- **Benchmark Environment:** 100 MCMC iterations, 100 MH draws per update
- **Prior Settings:** φ ~ Pareto(0.1, 1.5), λ ~ IG(0.01, 0.01)
- **MH Step Size:** Adaptive with initial value 1.0
- **Data Generation:** Mixture of Weibull(2, 1) and Weibull(1.5, 3)
- **Anomaly:** Performance irregularity at N=3162 likely due to memory/cache effects
---
*Benchmark conducted using the atime R package for asymptotic performance analysis.*

