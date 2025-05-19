# ClusterComponentUpdate Algorithm Analysis

## Algorithm Overview

This implements Neal's Algorithm 4 (for conjugate cases) and Algorithm 8 (for non-conjugate cases) for updating cluster assignments of data points in a Dirichlet Process Mixture Model (DPMM).

## Implementation Structure

The function is implemented with S3 dispatch based on the `dpObj` class:
- `ClusterComponentUpdate()`: Generic function.
- `ClusterComponentUpdate.conjugate()`: Handles DPMMs with conjugate base measures.
- `ClusterComponentUpdate.nonconjugate()`: Handles DPMMs with non-conjugate base measures.
- `ClusterComponentUpdate.hierarchical()`: Handles hierarchical Dirichlet Processes.

## Algorithmic Steps (Conjugate Case - `ClusterComponentUpdate.conjugate`)

Refers to Neal's Algorithm 4.
1. For each data point `i` from 1 to `n` (total number of data points):
   a. Temporarily remove data point `y[i,]` from its current cluster. Update `pointsPerCluster` for the current cluster label.
   b. Calculate probabilities for assigning `y[i,]` to each existing cluster `k`:
      `prob_existing[k] = pointsPerCluster[k] * Likelihood(mdObj, y[i,], clusterParams_k)`
      where `mdObj` is the mixing distribution object, and `clusterParams_k` are the parameters for cluster `k`.
   c. Calculate the probability for assigning `y[i,]` to a new cluster:
      `prob_new = alpha * predictiveArray[i]`
      where `alpha` is the concentration parameter and `predictiveArray[i]` is the predictive likelihood of `y[i,]` given the base measure.
   d. Normalize these probabilities (ensure they sum to 1, handle cases where all are zero).
   e. Sample a new cluster assignment (`newLabel`) for `y[i,]` based on these probabilities.
   f. Call `ClusterLabelChange()` to update `clusterLabels`, `pointsPerCluster`, `clusterParams`, and `numberClusters` based on `newLabel` and the `currentLabel`. This function handles creating a new cluster if `newLabel` indicates a new cluster, or removing the old cluster if it becomes empty.

## Algorithmic Steps (Non-conjugate Case - `ClusterComponentUpdate.nonconjugate`)

Refers to Neal's Algorithm 8.
1. For each data point `i` from 1 to `n`:
   a. Temporarily remove data point `y[i,]` from its current cluster. Update `pointsPerCluster`.
   b. Calculate probabilities for assigning `y[i,]` to each existing cluster `k`:
      `prob_existing[k] = pointsPerCluster[k] * Likelihood(mdObj, y[i,], clusterParams_k)`
   c. If the current cluster becomes empty after removing `y[i,]`, draw `m-1` auxiliary parameters from the prior `PriorDraw(mdObj, m-1)`. Combine these with the parameters of the (now empty) current cluster to form `m` sets of auxiliary parameters `aux`.
   d. If the current cluster is not empty, draw `m` auxiliary parameters `aux` from `PriorDraw(mdObj, m)`.
   e. Calculate probabilities for assigning `y[i,]` to new clusters using these auxiliary parameters:
      `prob_new_aux[j] = (alpha/m) * Likelihood(mdObj, y[i,], aux_j)` for `j` from 1 to `m`.
   f. Combine `prob_existing` and `prob_new_aux`. Normalize these probabilities.
   g. Sample a new cluster assignment (`newLabel`) for `y[i,]`.
   h. Call `ClusterLabelChange()` to update state. If `newLabel` corresponds to an auxiliary parameter, that auxiliary parameter becomes the parameter for the newly formed/assigned cluster.

## Algorithmic Steps (Hierarchical Case - `ClusterComponentUpdate.hierarchical`)

1. Iterates through each individual Dirichlet Process object (`indDP`) in the hierarchical structure.
2. For each `indDP`, calls the appropriate `ClusterComponentUpdate` method (conjugate or non-conjugate).
3. After updating components for an `indDP`, calls `DuplicateClusterRemove()` to merge any clusters that may have become identical during the update process.

## Key Data Structures Used

- `dpObj$data`: The input data (matrix).
- `dpObj$n`: Number of data points.
- `dpObj$alpha`: Concentration parameter of the DP.
- `dpObj$clusterLabels`: Vector storing cluster assignment for each data point.
- `dpObj$clusterParameters`: List of parameters for each cluster.
- `dpObj$numberClusters`: Current number of active clusters.
- `dpObj$pointsPerCluster`: Vector storing the count of data points in each cluster.
- `dpObj$mixingDistribution`: Object containing functions like `Likelihood()`, `PriorDraw()`, `Predictive()`.
- `dpObj$predictiveArray` (conjugate only): Stores precomputed predictive likelihoods.
- `dpObj$m` (non-conjugate only): Number of auxiliary parameters.

## Performance Considerations

1.  **Sequential Loop:** The primary loop iterates through each data point sequentially, making parallelization for this part challenging.
2.  **Likelihood Computations:** `Likelihood()` is called multiple times within the loop (for each existing cluster and for new/auxiliary clusters). Efficient likelihood calculation is crucial.
3.  **`ClusterLabelChange()` Overhead:** This function can involve re-allocating and resizing `clusterParameters` and `pointsPerCluster` if clusters are created or deleted, which can be costly.
4.  **`PriorDraw()` (non-conjugate):** Drawing `m` auxiliary parameters in each iteration adds computational load.
5.  **S3 Dispatch:** Minor overhead from S3 method dispatch, but likely less significant than the above points.

## Dependencies

- `Likelihood()` method for the specific mixing distribution.
- `Predictive()` method (conjugate case) for the mixing distribution.
- `PriorDraw()` method (non-conjugate case) for the mixing distribution.
- `ClusterLabelChange()` function (handles the logic of updating cluster assignments, parameters, and counts after a new label is chosen).
- `DuplicateClusterRemove()` function (hierarchical case).

