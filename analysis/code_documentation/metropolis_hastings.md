# Metropolis-Hastings Algorithm for Non-Conjugate Parameter Updates

When updating cluster parameters $\phi_k$ for a non-conjugate model, we sample from the posterior $p(\phi_k | D_k, G_0) \propto \left( \prod_{y_j \in D_k} F(y_j | \phi_k) \right) p_0(\phi_k)$, where $p_0(\phi_k)$ is the prior density from $G_0$.

Algorithm steps (for one cluster $k$ at MCMC iteration $t+1$):
1. Current state: $\phi_k^{(t)}$.
2. Propose a new state $\phi_k'$ from a proposal distribution $q(\phi_k' | \phi_k^{(t)})$. A common choice is a symmetric random walk, e.g., $\phi_k' \sim N(\phi_k^{(t)}, \Sigma_{prop})$.
3. Calculate the acceptance ratio $\mathcal{A}(\phi_k', \phi_k^{(t)})$:
   $$\mathcal{A}(\phi_k', \phi_k^{(t)}) = \min\left(1, \frac{p(\phi_k' | D_k, G_0) q(\phi_k^{(t)} | \phi_k')}{p(\phi_k^{(t)} | D_k, G_0) q(\phi_k' | \phi_k^{(t)})}\right)$$
   If $q$ is symmetric (e.g., Gaussian random walk), $q(\phi_k^{(t)} | \phi_k') = q(\phi_k' | \phi_k^{(t)})$, so these terms cancel:
   $$\mathcal{A}(\phi_k', \phi_k^{(t)}) = \min\left(1, \frac{\left( \prod_{y_j \in D_k} F(y_j | \phi_k') \right) p_0(\phi_k')}{\left( \prod_{y_j \in D_k} F(y_j | \phi_k^{(t)}) \right) p_0(\phi_k^{(t)})}\right)$$
4. Draw $u \sim U(0,1)$.
5. If $u < \mathcal{A}(\phi_k', \phi_k^{(t)})$, set $\phi_k^{(t+1)} = \phi_k'$ (accept).
6. Else, set $\phi_k^{(t+1)} = \phi_k^{(t)}$ (reject).

Implementation typically uses log-probabilities for numerical stability.
