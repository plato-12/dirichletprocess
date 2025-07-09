# Posterior Draw: Beta Parameters (Non-Conjugate via Metropolis-Hastings)

For a Beta likelihood $F(y | \mu_y, \tau, M_{max})$ with non-conjugate priors on $(\mu_y, \tau)$ (e.g., $U(0,M_{max})$ for $\mu_y$, $\text{Gamma}(a_\tau, b_\tau)$ for $\tau$), we use Metropolis-Hastings to sample from the posterior.

Given data $D_k = \{y_j : c_j=k\}$ for cluster $k$.
Target posterior density:
$$p(\mu_y, \tau | D_k, G_0) \propto \left( \prod_{y_j \in D_k} F(y_j | \mu_y, \tau, M_{max}) \right) p_0(\mu_y) p_0(\tau)$$
where $p_0(\cdot)$ are the prior densities from $G_0$.

Metropolis-Hastings steps for $(\mu_y, \tau)$:
1. Initialize $(\mu_y^{(t)}, \tau^{(t)})$.
2. Propose $(\mu_y', \tau')$ from $q((\mu_y', \tau') | (\mu_y^{(t)}, \tau^{(t)}))$. E.g., independent random walks for $\mu_y$ (clipped to $[0, M_{max}]$) and $\log(\tau)$ or $\tau$ (ensuring positivity).
3. Calculate acceptance ratio $\mathcal{A}$ (see `metropolis_hastings.md` for general formula).
4. Accept or reject the proposal to get $(\mu_y^{(t+1)}, \tau^{(t+1)})$.

The `dirichletprocess` package uses specific proposal mechanisms (`MhParameterProposal.beta`) and prior densities (`PriorDensity.beta`) within the MH algorithm.
