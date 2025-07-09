# ClusterLabelChange() Function Analysis

## Algorithm Overview

The `ClusterLabelChange()` S3 generic method is a utility function called within `ClusterComponentUpdate()`. After a new cluster label (`newLabel`) has been sampled for a data point `i` (which was previously in `currentLabel`), this function updates the state of the Dirichlet Process object (`dpObj`) to reflect this change. This includes updating cluster assignments, counts of points per cluster, the cluster parameters list, and the total number of clusters.

## Implementation Structure

S3 dispatch is used based on the `dpObj`'s conjugacy status:
- `ClusterLabelChange.conjugate()`
- `ClusterLabelChange.nonconjugate()`
The two methods handle slightly different logic regarding how new cluster parameters are incorporated, especially when a new cluster is formed from an auxiliary parameter in the non-conjugate case.

## Algorithmic Steps (Common Logic and `ClusterLabelChange.conjugate`)

Inputs: `dpObj`, data point index `i`, `newLabel`, `currentLabel`, `aux` (auxiliary parameters, primarily for non-conjugate).
1. Retrieve current state: `pointsPerCluster`, `clusterLabels`, `clusterParameters`, `numLabels`, `mdObj` from `dpObj`.
2. **Case 1: `newLabel` corresponds to an existing cluster (`newLabel <= numLabels`)**
   a. Increment count for the new cluster: `pointsPerCluster[newLabel] = pointsPerCluster[newLabel] + 1`.
   b. Update the data point's label: `clusterLabels[i] = newLabel`.
   c. Check if the `currentLabel` cluster has become empty (`pointsPerCluster[currentLabel] == 0`):
      i. If empty: Decrement `numLabels`. Remove the `currentLabel` entry from `pointsPerCluster` and from each parameter array in `clusterParameters`.
      ii. Re-index `clusterLabels`: Any label greater than `currentLabel` is decremented by 1 to fill the gap.
3. **Case 2: `newLabel` corresponds to a new cluster being formed (`newLabel > numLabels`)**
   a. **If `currentLabel` became empty (`pointsPerCluster[currentLabel] == 0`)** (This means the point was the last one in its old cluster, and it's forming a new cluster by itself):
      i.  **Conjugate**: Draw new parameters for this now re-purposed cluster slot: `post_draw = PosteriorDraw(mdObj, dpObj$data[i, , drop = FALSE])`. Update `clusterParameters[[...]][,,currentLabel]` with `post_draw`.
      ii. **Non-conjugate**: Assign parameters from the chosen auxiliary parameter: `clusterParameters[[...]][,,currentLabel] = aux[[...]][,, newLabel - numLabels]` (where `newLabel - numLabels` is the index into the `aux` array).
      iii.Increment its count: `pointsPerCluster[currentLabel] = pointsPerCluster[currentLabel] + 1` (it becomes 1).
      iv. The data point `i` is already implicitly assigned to `currentLabel` if this path is taken (its label doesn't need to change if it's re-using the slot of its just-emptied cluster for the new parameters).
   b. **Else (`currentLabel` did not become empty OR it was already a new cluster)**:
      i.  Assign `clusterLabels[i] = numLabels + 1` (or a re-indexed version if an empty slot was used).
      ii. Increment `numLabels`.
      iii.Append 1 to `pointsPerCluster`.
      iv. **Conjugate**: Draw parameters for the new cluster: `post_draw = PosteriorDraw(mdObj, dpObj$data[i, , drop = FALSE])`. Append `post_draw` to `clusterParameters`.
      v.  **Non-conjugate**: Take parameters from the chosen auxiliary parameter `aux[[...]][,, newLabel - numLabels]` and append to `clusterParameters`.
4. Update `dpObj` with the modified `pointsPerCluster`, `clusterLabels`, `clusterParameters`, and `numLabels`.

## `ClusterLabelChange.nonconjugate` Specifics

- The main difference lies in step 3.b.iv/v and 3.a.ii. When a new cluster is created (either by splitting off or by a point moving to an auxiliary parameter slot), its parameters are taken from the `aux` list (the auxiliary parameters drawn in `ClusterComponentUpdate.nonconjugate`). The index `newLabel - numLabels` is used to select which set of auxiliary parameters from `aux` to use.

## Key Data Structures Manipulated

- `dpObj$pointsPerCluster`: Counts are incremented/decremented. May be resized if a cluster is removed or added.
- `dpObj$clusterLabels`: Assignment for data point `i` is changed.
- `dpObj$clusterParameters`: Parameters for an existing cluster might be removed (if it becomes empty and its slot is reused or removed). New parameters are added if a new cluster is genuinely created (not reusing an empty slot).
- `dpObj$numberClusters`: Incremented or decremented.
- `aux` (Input, non-conjugate only): List of auxiliary parameter sets.

## Performance Considerations

1.  **Array/List Manipulation**: The primary costs come from potentially resizing/modifying `clusterParameters` (a list of arrays) and `pointsPerCluster` (a vector). Appending or removing elements from lists/vectors in R can involve copying data.
2.  **Re-indexing Labels**: If a cluster is removed, `clusterLabels` needs to be re-indexed, which is a vector operation over `n` elements.
3.  **`PosteriorDraw` Call (Conjugate, new cluster from empty old one)**: In the conjugate case, if a point leaves its cluster making it empty, and then forms a *new* cluster *by itself reusing that slot*, `PosteriorDraw` is called for that single data point. This is generally fast.
4.  **Parameter Assignment (Non-conjugate)**: Simply assigns parameters from the `aux` list, which is efficient.
5.  **Frequency of Calls**: Called once for every data point in every MCMC iteration within `ClusterComponentUpdate`.

## Dependencies

- `PosteriorDraw()` (conjugate case, for initializing a new cluster that reuses an emptied slot, or for initializing a completely new cluster slot).
- (Implicitly) `PriorDraw()` for generating `aux` parameters in the calling function (`ClusterComponentUpdate.nonconjugate`).

