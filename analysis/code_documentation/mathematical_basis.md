# Mathematical Foundation of dirichletprocess Algorithms

## Dirichlet Process Mixture Model (DPMM)

The model assumes data $y_i$ are generated from a mixture model where the mixing distribution $G$ itself is drawn from a Dirichlet Process (DP). The hierarchical representation is:

$$y_i | \theta_i \sim F(\cdot | \theta_i)$$
$$\theta_i | G \sim G$$
$$G | G_0, \alpha \sim DP(G_0, \alpha)$$

Where:
- $y_i$ is the $i$-th observation.
- $F(\cdot | \theta_i)$ is the likelihood function for an observation given parameters $\theta_i$ (e.g., Gaussian $N(y_i | \mu_i, \sigma_i^2)$, Beta $Beta(y_i | a_i, b_i)$).
- $\theta_i$ are the parameters specific to the $i$-th observation.
- $G$ is a random probability measure (the mixing distribution) drawn from a Dirichlet Process.
- $G_0$ is the base measure, representing the prior belief about the distribution of parameters $\theta$. It is the expected value of $G$, $E[G] = G_0$.
- $\alpha > 0$ is the concentration parameter. Larger $\alpha$ leads to $G$ being closer to $G_0$ and implies a larger number of distinct parameter values (clusters) a priori.

Due to the discreteness of draws from a DP, many $\theta_i$ will share common values. Let $\{\phi_k\}_{k=1}^K$ be the set of unique parameter values, where $K$ is the number of clusters.

## Chinese Restaurant Process (CRP) / Gibbs Sampling for Cluster Assignments (Neal's Algorithm 3)

The Chinese Restaurant Process (CRP) formulation provides a constructive way to understand the clustering property of the DP. For Gibbs sampling, we integrate out $G$. The conditional posterior probability for assigning data point $y_i$ to cluster $k$ (parameterized by $\phi_k$), given all other assignments $c_{-i}$ and current cluster parameters, is (Neal, 2000, Algorithm 3):

$$P(c_i = k | c_{-i}, \mathbf{y}, \alpha, G_0, \{\phi_j\}_{j \neq k}) \propto \begin{cases}
n_{k,-i} \cdot F(y_i | \phi_k) & \text{for an existing cluster } k \text{ (where } n_{k,-i} > 0) \\
\alpha \cdot \int F(y_i | \phi) dG_0(\phi) & \text{for a new cluster (i.e., } k = K_{new})
\end{cases}$$

Where:
- $c_i$ is the cluster assignment for data point $y_i$.
- $c_{-i}$ denotes the cluster assignments for all data points except $y_i$.
- $n_{k,-i}$ is the number of data points in cluster $k$, excluding point $y_i$.
- $F(y_i | \phi_k)$ is the likelihood of data point $y_i$ given the parameters $\phi_k$ of cluster $k$.
- $\alpha$ is the concentration parameter.
- $\int F(y_i | \phi) dG_0(\phi)$ is the marginal likelihood of $y_i$ under the base measure $G_0$, often called the predictive likelihood. This is used when proposing a new cluster, whose parameters $\phi_{K_{new}}$ would be drawn from $G_0$ (or a posterior based on $y_i$ if conjugate).

This is Neal's (2000) Algorithm 3. Algorithm 2 is similar but marginalizes over $\phi_k$ as well, which is common for fully conjugate models.

## Cluster Parameter Updates (Neal's Algorithm 3)

After re-assigning all $c_i$, the parameters $\phi_k$ for each cluster $k$ (now containing points $D_k = \{y_i : c_i = k\}$) are updated.

### Conjugate Case

If $G_0$ is conjugate to the likelihood $F$, the parameters $\phi_k$ for each cluster $k$ are drawn from their posterior distribution:

$$\phi_k | D_k, G_0 \sim p(\phi_k | D_k, G_0)$$

This posterior is analytically known. For example, if $F$ is Gaussian and $G_0$ is Normal-Inverse-Gamma, the posterior $p(\phi_k | D_k, G_0)$ is also Normal-Inverse-Gamma with updated hyperparameters.

### Non-conjugate Case (Neal's Algorithm 8 for assignments, Metropolis-Hastings for parameters)

When $G_0$ is not conjugate to $F$, sampling $\phi_k$ directly from the posterior is not possible. Instead, MCMC methods like Metropolis-Hastings are used. The target density is:

$$p(\phi_k | D_k, G_0) \propto \left( \prod_{y_j \in D_k} F(y_j | \phi_k) \right) p_0(\phi_k)$$

where $p_0(\phi_k)$ is the prior density of $\phi_k$ from $G_0$.
For cluster assignments in non-conjugate cases, Neal's Algorithm 8 is often used, which introduces $m$ auxiliary parameters drawn from $G_0$ to help propose new clusters, avoiding the direct calculation of $\int F(y_i | \phi) dG_0(\phi)$.

## Concentration Parameter $\alpha$ Update (Escobar & West, 1995)

The concentration parameter $\alpha$ can be given a prior (e.g., Gamma distribution $p(\alpha) = \text{Gamma}(\alpha | a_\alpha, b_\alpha)$) and updated using an auxiliary variable method (Escobar and West, 1995). Given $K$ current clusters and $N$ data points, the update involves drawing an auxiliary variable $\eta \sim \text{Beta}(\alpha+1, N)$. The posterior for $\alpha$ is then a mixture of two Gamma distributions:

$$p(\alpha | K, N, \eta, \text{priors}) = \pi_\eta \text{Gamma}(\alpha | a_\alpha+K, b_\alpha-\log(\eta)) + (1-\pi_\eta) \text{Gamma}(\alpha | a_\alpha+K-1, b_\alpha-\log(\eta))$$

Where the mixture weight $\pi_\eta$ depends on $a_\alpha, K, N, b_\alpha,$ and $\log(\eta)$.

