# ClusterParameterUpdate Algorithm Analysis

## Algorithm Overview

This function updates the parameters for each existing cluster based on the data points currently assigned to it. It's a crucial step in the MCMC fitting process of a Dirichlet Process Mixture Model (DPMM).

## Implementation Structure

The function is implemented with S3 dispatch based on the `dpObj` class (specifically, its conjugacy status):
- `ClusterParameterUpdate()`: Generic function.
- `ClusterParameterUpdate.conjugate()`: Handles DPMMs with conjugate base measures.
- `ClusterParameterUpdate.nonconjugate()`: Handles DPMMs with non-conjugate base measures.

Note: There isn't a specific `ClusterParameterUpdate.hierarchical()` method explicitly shown in `cluster_parameter_update.R`. Hierarchical parameter updates are likely handled within `GlobalParameterUpdate` and updates to individual DPs within the hierarchy would use their respective conjugate/non-conjugate methods. The proposal mentions `GlobalParameterUpdate` for hierarchical models (Section 8.1 of priyanshuTiwari_dirichletprocess.pdf).

## Algorithmic Steps (Conjugate Case - `ClusterParameterUpdate.conjugate`)

1. For each existing cluster `i` from 1 to `numLabels` (current number of clusters):
   a. Identify all data points `pts` assigned to cluster `i` (i.e., where `clusterLabels == i`).
   b. Draw new parameters for cluster `i` from the posterior distribution: `post_draw = PosteriorDraw(mdobj, pts)`.
      `mdobj` is the mixing distribution object. `PosteriorDraw` leverages conjugacy to sample directly from the analytical posterior given the data `pts` assigned to the cluster.
   c. Update `clusterParams` for cluster `i` with these `post_draw` values.

## Algorithmic Steps (Non-conjugate Case - `ClusterParameterUpdate.nonconjugate`)

1. For each existing cluster `i` from 1 to `numLabels`:
   a. Identify all data points `pts` assigned to cluster `i`.
   b. Get the current parameters for cluster `i` to use as a starting position for the Metropolis-Hastings sampler: `start_pos`.
   c. Draw `mhDraws` samples from the posterior distribution using Metropolis-Hastings: `parameter_samples = PosteriorDraw(mdobj, pts, mhDraws, start_pos = start_pos)`.
      `PosteriorDraw` for non-conjugate cases internally uses `MetropolisHastings()` which involves proposing new parameters (`MhParameterProposal()`), calculating likelihoods (`Likelihood()`), prior densities (`PriorDensity()`), and acceptance probabilities.
   d. Update `clusterParams` for cluster `i` with the last sample from `parameter_samples`.
   e. (Optionally) Calculate and store the acceptance ratio for the MH step for diagnostics.

The internal function `cluster_parameter_update` (not exported, and seemingly a conceptual or alternative implementation) also outlines a similar process where it iterates unique clusters and calls `PosteriorDraw` for data belonging to each.

## Key Data Structures Used

- `dpObj$data`: The input data.
- `dpObj$clusterLabels`: Vector of cluster assignments.
- `dpObj$clusterParameters`: List of parameters for each cluster (this is what gets updated).
- `dpObj$numberClusters`: Current number of active clusters.
- `dpObj$mixingDistribution`: Object containing `PosteriorDraw()`, `PriorDraw()`, `Likelihood()`, `PriorDensity()`, `MhParameterProposal()` methods.
- `dpObj$mhDraws` (non-conjugate only): Number of Metropolis-Hastings iterations for parameter sampling.

## Performance Considerations

1.  **Looping over Clusters:** The algorithm iterates through each active cluster.
2.  **Data Subsetting:** For each cluster, it subsets the data points belonging to it.
3.  **`PosteriorDraw()` Cost:**
    -   **Conjugate:** Relatively efficient as it involves direct sampling from a known posterior distribution. The complexity depends on calculating posterior parameters from data in the cluster.
    -   **Non-conjugate:** Can be significantly more expensive due to the iterative nature of `MetropolisHastings()`. This involves `mhDraws` iterations of proposing parameters, and repeatedly calculating likelihoods and prior densities.
4.  **Likelihood and Prior Density Calls (Non-conjugate):** Inside `MetropolisHastings`, `Likelihood()` and `PriorDensity()` are called for each of the `mhDraws` iterations per cluster.

## Dependencies

- `PosteriorDraw()` method for the specific mixing distribution (which itself depends on other methods for non-conjugate cases, such as `MetropolisHastings`, `MhParameterProposal`, `Likelihood`, `PriorDensity`).

