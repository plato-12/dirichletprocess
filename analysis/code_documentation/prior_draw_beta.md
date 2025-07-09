# Prior Draw: Beta Parameters from $G_0$

For a Beta likelihood $F(y | \mu_y, \tau, M_{max})$ (mean $\mu_y \in [0, M_{max}]$, precision $\tau$), a typical non-conjugate $G_0$ might define independent priors for $\mu_y$ and $\tau$.

Hyperparameters for $G_0$ could be:
- For $\mu_y$: Often $U(0, M_{max})$.
- For $\tau$: e.g., $\text{Gamma}(a_\tau, b_\tau)$ or $\text{Inv-Gamma}(a_\tau, b_\tau)$. The `dirichletprocess` R code for `beta_uniform_gamma` uses an Inverse Gamma prior for `nu` (scale parameter, $\tau = \text{nu}$ in that context), meaning `nu` ~ InvGamma(alpha0, beta0).

Sampling procedure (example using $U(0,M_{max})$ for $\mu_y$ and $\text{Gamma}(a_\tau,b_\tau)$ for $\tau$ as precision):
1. Draw mean $\mu_y$:
   $$\mu_y \sim U(0, M_{max})$$
2. Draw precision $\tau$:
   $$\tau \sim \text{Gamma}(a_\tau, b_\tau)$$

The pair $(\mu_y, \tau)$ is one draw from $G_0$.
