# Hierarchical Dirichlet Process (HDP) - Mathematical Details

The HDP (Teh et al., 2005) allows sharing of statistical strength across multiple related groups, where each group $j$ has its own DP mixture model, but these DPs share a common, global base measure $G_0$, which is itself drawn from a DP.

## Model Specification

For group $j=1, \dots, J$:
$$y_{ji} | \theta_{ji} \sim F(\cdot | \theta_{ji}) \quad \text{(Observation } i \text{ in group } j)$$
$$\theta_{ji} | G_j \sim G_j \quad \text{(Parameters for observation } ji)$$
$$G_j | G_0, \alpha_0 \sim DP(G_0, \alpha_0) \quad \text{(Group-specific DP, } \alpha_0 \text{ is concentration for } G_j)$$
$$G_0 | H, \gamma \sim DP(H, \gamma) \quad \text{(Global DP, } \gamma \text{ is top-level concentration})$$

Where:
- $F$ is the likelihood function.
- $G_j$ is the random measure for group $j$.
- $G_0$ is the global random measure, serving as the base measure for each $G_j$. It is discrete: $G_0 = \sum_{k=1}^{\infty} \beta_k \delta_{\phi_k}$.
- $\phi_k \sim H$ are the global atoms (parameters) drawn from a hyperprior $H$.
- $\boldsymbol{\beta} = (\beta_1, \beta_2, \dots)$ are the global stick-breaking weights from $DP(H, \gamma)$, i.e., $\beta_k = v_k \prod_{l<k}(1-v_l)$ with $v_k \sim \text{Beta}(1, \gamma)$.
- $\alpha_0$ is the concentration parameter for group-level DPs (can be different per group, $\alpha_j$).

## Key Update Steps in HDP Gibbs Sampler (Simplified, based on Neal's methods where applicable)

### 1. Update Data Point Assignments $c_{ji}$ (within each group $j$)
For each data point $y_{ji}$ in group $j$, assign it to a local cluster $l$ (with params $\theta_{jl}$) or a new local cluster using the CRP specific to $G_j$. The probability of assigning to an existing local cluster $l$ (with $n_{jl,-i}$ points) is proportional to $n_{jl,-i} F(y_{ji} | \theta_{jl})$. The probability of assigning to a new local cluster is proportional to $\alpha_0 \int F(y_{ji} | \theta) dG_0(\theta)$. The parameters $\theta_{jl}$ are themselves draws from $G_0$.

### 2. Update Local Cluster Parameters $\theta_{jl}$ (within each group $j$)
For each local cluster $l$ in group $j$, its parameters $\theta_{jl}$ are effectively a draw from the global $G_0$. In practice, these are often linked to one of the global atoms $\phi_k$. The update involves re-sampling which global atom $\phi_k$ this local cluster type $\theta_{jl}$ corresponds to, with probability proportional to $\beta_k \prod_{y \in D_{jl}} F(y | \phi_k)$.

### 3. Update Global Atoms $\phi_k$
Each global atom $\phi_k$ is updated based on all data points $D_k = \{ y_{ji} : \theta_{ji} \text{ is an instance of } \phi_k \}$ that are associated with it, across all groups. The update is from the posterior:
$$\phi_k | D_k, H \sim p(\phi_k | D_k, H) \propto \left( \prod_{y \in D_k} F(y | \phi_k) \right) p_H(\phi_k)$$
where $p_H(\phi_k)$ is the prior density from $H$. This step uses standard conjugate updates or MH if $H$ is not conjugate to $F$.

### 4. Update Global Stick-Breaking Weights $\boldsymbol{\beta}$ (and potentially add/remove atoms from $G_0$)
Let $m_k$ be the number of distinct local cluster types (tables) across all groups that are instances of the global atom $\phi_k$. Let $K_0$ be the number of currently instantiated global atoms.
The posterior for the weights $(\beta_1, \dots, \beta_{K_0}, \beta_{new_mass})$ is Dirichlet:
$$(\beta_1, \dots, \beta_{K_0}, \sum_{l=K_0+1}^{\infty} \beta_l) \sim \text{Dirichlet}(m_1, \dots, m_{K_0}, \gamma)$$
The sum $\sum_{l=K_0+1}^{\infty} \beta_l$ is the mass for new global atoms. New atoms $\phi_{K_0+1}, \dots$ can be drawn from $H$, and their weights derived by further stick-breaking on this remaining mass, often with truncation.

### 5. Update Concentration Parameters $\alpha_0$ (group-level) and $\gamma$ (global-level)
- Each $\alpha_j$ (or common $\alpha_0$) is updated using the Escobar & West method, based on the number of local clusters $K_j$ and data points $N_j$ within group $j$.
- The global $\gamma$ is updated using Escobar & West, based on the number of distinct global atoms $K_0$ currently used by any table, and the total number of tables $M = \sum_k m_k$ across all groups.

This is a simplified overview. Specific algorithms (e.g., Neal's Algorithm 8 for HDPs, or direct assignment methods) have more detailed steps for managing assignments to global atoms and sampling new atoms.
