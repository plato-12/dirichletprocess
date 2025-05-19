# ClusterComponentUpdate Algorithm Analysis

## Algorithm Overview

This implements Neal's Algorithm 4 (for conjugate cases) and Algorithm 8 (for non-conjugate cases).

## Implementation Structure

The function is implemented with S3 dispatch:
- `ClusterComponentUpdate()`: Generic function
- `ClusterComponentUpdate.conjugate()`: For conjugate mixtures
- `ClusterComponentUpdate.nonconjugate()`: For non-conjugate mixtures
- `ClusterComponentUpdate.hierarchical()`: For hierarchical models

## Algorithmic Steps (Conjugate Case)

1. For each data point i:
   a. Remove point from current cluster
   b. Calculate probabilities for each existing cluster
   c. Calculate probability for a new cluster
   d. Sample new cluster assignment
   e. Update cluster assignments

## Algorithmic Steps (Non-conjugate Case)

1. For each data point i:
   a. Remove point from current cluster
   b. Calculate probabilities for each existing cluster
   c. Draw auxiliary parameters
   d. Calculate probabilities for new clusters with auxiliary parameters
   e. Sample new cluster assignment
   f. Update cluster assignments

## Key Data Structures

- `clusterLabels`: Vector of cluster assignments
- `pointsPerCluster`: Count of points in each cluster
- `clusterParams`: Parameters for each cluster

## Performance Considerations

1. The algorithm is inherently sequential (loop through each data point)
2. Likelihood calculations are performed repeatedly
3. Memory allocations occur during cluster creation/deletion

