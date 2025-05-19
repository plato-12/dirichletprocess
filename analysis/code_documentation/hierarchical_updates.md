# Hierarchical Dirichlet Process (HDP) Specific Updates Analysis

## Algorithm Overview

For Hierarchical Dirichlet Process models, in addition to updating parameters within each individual DP (`indDP`), global parameters that link these DPs must also be updated. This includes the global base measure $G_0$ (represented by its atoms and stick-breaking weights), the parameters of these global atoms, and the global concentration parameter $\gamma$.

Key functions involved:
- `GlobalParameterUpdate()`: Updates the parameters $\phi_k$ of the atoms in the global base measure $G_0$.
- `UpdateG0()`: Updates the discrete global base measure $G_0$ itself by resampling its stick-breaking weights and potentially atoms. This also updates the mixing distributions of individual DPs (`indDP`) to draw from this new $G_0$.
- `UpdateGamma()` (via `update_concentration`): Updates the top-level concentration parameter $\gamma$ for $G_0$.

## `GlobalParameterUpdate.hierarchical()` (global_parameter_update.R)

### Algorithm Overview
This function updates the parameters $\phi_k$ of the shared atoms in the global base measure $G_0$. Each $\phi_k$ is updated based on all data points from all individual DPs that are currently associated with that specific global atom.

### Algorithmic Steps
1. Identify unique global atoms currently in use: `global_labels = unique(unlist(lapply(indDPs, function(dp) match(dp$clusterParameters, G0_atoms))))` (conceptual representation).
2. For each unique global atom $\phi_k$ (indexed by `global_labels[i]`):
   a. Collect all data points `total_pts` from all individual DPs that are currently assigned to local clusters whose parameters are instances of this global atom $\phi_k$.
      - This involves iterating through each `indDP`.
      - Identifying which local clusters in `indDP` correspond to the current `global_labels[i]`.
      - Subsetting the data from `indDP` belonging to these local clusters.
   b. Draw a new parameter value for this global atom from its posterior distribution given `total_pts` and the hyperprior $H$ (the prior for $G_0$'s atoms): `new_param = PosteriorDraw(mixingDistribution_of_G0, total_pts, ...)`.
      - The `mixingDistribution` of $G_0$ defines the kernel and prior for $\phi_k$.
      - This step uses `PosteriorDraw`, which will be conjugate or non-conjugate based on $H$.
   c. Update the `globalParameters` list (which stores the $\phi_k$ values) with `new_param`.
   d. **Crucially**: Update the parameters of all local clusters across all `indDPs` that were instances of the *old* $\phi_k$ to now be instances of the *new* (updated) $\phi_k$.
3. Update the `mixingDistribution$theta_k` within each `indDP` to reflect the (potentially changed) `globalParameters`.

## `UpdateG0()` (update_g0.R)

### Algorithm Overview
This function updates the discrete global base measure $G_0 = \sum \beta_k \delta_{\phi_k}$ by resampling the global stick-breaking weights $\beta_k$ and using the current set of global atoms $\phi_k$ (which might have been updated by `GlobalParameterUpdate`). It then updates the individual DPs to draw their cluster parameters from this new $G_0$.

### Algorithmic Steps
1. Get current `globalParameters` ($\{\phi_k\}$) and global concentration `gamma`.
2. Determine usage of global atoms: Count how many distinct local cluster types (tables) across all `indDPs` are associated with each global atom $\phi_k$. This forms `globalParamTable$Freq`.
3. Resample global stick-breaking weights $\beta_k$:
   a. Use a Dirichlet distribution draw: `dirichlet_draws = gtools::rdirichlet(1, c(globalParamTable$Freq, dpobjlist$gamma))`.
      The counts `globalParamTable$Freq` are the counts of tables (unique cluster parameter values across all groups) that map to each global atom. `dpobjlist$gamma` is the concentration parameter for the top-level DP generating $G_0$.
   b. Truncate the stick-breaking process for the 'new' part of $G_0$: Generate `numBreaks` new sticks using `StickBreaking(dpobjlist$gamma + numTables, numBreaks)`, where `numTables` is sum of `globalParamTable$Freq`.
   c. Combine: `sticks = c(dirichlet_draws[1:numGlobalAtomsUsed], dirichlet_draws[numGlobalAtomsUsed+1] * new_sticks_from_StickBreaking)`.
4. Generate new atoms if needed for the newly created sticks from the truncation:
   `priorDraws = PriorDraw(mixingDistribution_of_G0, numBreaks)`.
5. Combine existing global atoms with newly drawn ones: `postParams` becomes `c(current_global_atoms, priorDraws)`.
6. For each individual DP (`indDP[[i]]`):
   a. Generate new local stick-breaking weights `pi_k^(j)` for this DP using the new global sticks `beta_k` (now called `sticks`) and the `indDP[[i]]$alpha`: `newGJ = draw_gj(indDP[[i]]$alpha, sticks)`.
   b. Update `indDP[[i]]$mixingDistribution$pi_k = newGJ`.
   c. Update `indDP[[i]]$mixingDistribution$theta_k = postParams` (all indDPs now point to the same updated set of global atoms and their parameters).
7. Update `dpobjlist$globalParameters = postParams` and `dpobjlist$globalStick = sticks`.

## `UpdateGamma()` (called via `UpdateAlpha.hierarchical` -> `update_concentration` for the global gamma)

### Algorithm Overview
Updates the top-level concentration parameter $\gamma$ for the DP that generates $G_0$, using the same Escobar & West method as `UpdateAlpha.default`, but applied to the 'tables' and 'customers' at the global level.

### Algorithmic Steps
1. Count `numParams`: Number of unique global atoms currently used (i.e., number of distinct $\phi_k$ values that have at least one local cluster mapped to them).
2. Count `numTables`: Total number of distinct local cluster types across all `indDPs` (i.e., sum of the number of unique cluster parameters in each group, before collapsing to global atoms, effectively this is the sum of `globalParamTable$Freq` from `UpdateG0`).
3. Call `update_concentration(current_gamma, numTables, numParams, gammaPriors)`.

## Key Data Structures Used

- `dpobjlist$indDP`: List of individual Dirichlet Process objects.
- `dpobjlist$globalParameters`: List of arrays storing the parameters $\{\phi_k\}$ of the global atoms.
- `dpobjlist$globalStick`: Vector of stick-breaking weights $\{\beta_k\}$ for $G_0$.
- `dpobjlist$gamma`: Global concentration parameter.
- `dpobjlist$gammaPriors`: Prior for `gamma`.
- `indDP[[i]]$clusterParameters`: Parameters of local clusters.
- `indDP[[i]]$clusterLabels`: Local cluster assignments.
- `indDP[[i]]$mixingDistribution$pi_k`: Local stick-breaking weights (pointing to global atoms).
- `indDP[[i]]$mixingDistribution$theta_k`: Reference to global atoms $\{\phi_k\}$.

## Performance Considerations

1.  **`GlobalParameterUpdate`**:
    -   Iterates over unique global atoms.
    -   For each global atom, iterates over all `indDPs` and their local clusters to gather data.
    -   `PosteriorDraw()` for each global atom can be expensive, especially if non-conjugate and data pooled from many groups is large.
2.  **`UpdateG0`**:
    -   Counting global atom usage involves iterating through `indDPs` and their cluster parameters.
    -   `gtools::rdirichlet` and `StickBreaking` are generally efficient.
    -   `PriorDraw` for new global atoms depends on `numBreaks`.
    -   The loop to update `pi_k` for each `indDP` involves `draw_gj`, which itself has a loop.
3.  **Data Aggregation**: Collecting all data points associated with a global atom (`total_pts` in `GlobalParameterUpdate`) can be memory and time-consuming if not done efficiently.

## Dependencies

- All core DP functions for individual DPs (`PosteriorDraw`, `PriorDraw`, `Likelihood` for the kernel of $\phi_k$).
- `gtools::rdirichlet`.
- `StickBreaking()`, `draw_gj()`.
- `update_concentration()` (for `UpdateGamma`).
- `true_cluster_labels()` for correctly matching multi-dimensional parameters.

