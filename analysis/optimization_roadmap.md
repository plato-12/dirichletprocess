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
