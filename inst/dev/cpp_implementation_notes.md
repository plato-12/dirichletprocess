# C++ Implementation Notes for dirichletprocess

## Overview

This document describes the C++ implementations of core sampling algorithms for the dirichletprocess package.

## Implemented Functions

### Normal Distribution (Conjugate Case)

1. **PriorDraw**: `normal_prior_draw_cpp(priorParams, n)`
   - Draws from the Normal-Inverse-Gamma prior distribution
   - Parameters: mu0, kappa0, alpha0, beta0
   - Returns: List with mu and sigma arrays

2. **PosteriorDraw**: `normal_posterior_draw_cpp(priorParams, x, n)`
   - Draws from the posterior distribution given data
   - Computes posterior parameters first, then samples
   - Returns: List with mu and sigma arrays

3. **ClusterComponentUpdate**: `conjugate_cluster_component_update_cpp(dpObj)`
   - Implements Chinese Restaurant Process (Neal's Algorithm 4)
   - Updates cluster assignments for all data points
   - Returns: Updated cluster labels, sizes, parameters, and count

4. **ClusterParameterUpdate**: `conjugate_cluster_parameter_update_cpp(dpObj)`
   - Updates parameters for each cluster using posterior draws
   - Returns: Updated cluster parameters

## Key Implementation Details

### Array Indexing
- **Important**: C++ uses 0-based indexing while R uses 1-based indexing
- Cluster labels must be converted: `dp$clusterLabels - 1` before passing to C++
- Convert back after: `result$clusterLabels + 1`

### Random Number Generation
- Uses R's RNG through Rcpp (R::rnorm, R::rgamma, etc.)
- Set seed in R before calling C++ functions for reproducibility

### Memory Management
- Leverages Rcpp's automatic memory management
- No manual memory allocation/deallocation needed in most cases

## Usage Example

```r
# Load the package
library(dirichletprocess)

# Generate data
set.seed(123)
y <- c(rnorm(30, -3, 0.5), rnorm(30, 3, 0.5))

# Initialize DP
dp <- DirichletProcessGaussian(y)

# Convert to 0-indexed for C++
dp$clusterLabels <- dp$clusterLabels - 1

# Run one iteration of the sampler
# Step 1: Update cluster assignments
result <- conjugate_cluster_component_update_cpp(dp)
dp$clusterLabels <- result$clusterLabels
dp$pointsPerCluster <- result$pointsPerCluster
dp$numberClusters <- result$numberClusters
dp$clusterParameters <- result$clusterParameters

# Step 2: Update cluster parameters
dp$clusterParameters <- conjugate_cluster_parameter_update_cpp(dp)

# Convert back to 1-indexed
dp$clusterLabels <- dp$clusterLabels + 1
```

## Performance

Initial benchmarks show:
- Likelihood calculations: ~1.6x faster than R
- Full sampling iterations: Expected 10-100x speedup for large datasets

## Testing

Run tests with:
```r
devtools::test(filter = "normal-cpp")
```

## Future Extensions

The architecture supports adding:
- Beta distribution (non-conjugate)
- Weibull distribution (non-conjugate)
- Multivariate Normal distribution
- Hierarchical models

Each would follow the same pattern:
1. Implement distribution class inheriting from `MixingDistribution`
2. Implement DP class for conjugate/non-conjugate case
3. Export functions via Rcpp
4. Add comprehensive tests
