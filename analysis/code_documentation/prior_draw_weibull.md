# Prior Draw: Weibull Parameters ($k, \lambda$) from $G_0$

For a Weibull likelihood $F(y | k, \lambda)$ (shape $k$, scale $\lambda$), $G_0$ might define independent priors for $k$ and $\lambda$.

Hyperparameters for $G_0$ (example from `weibull_uniform_gamma.R`):
- For shape $k$: $U(0, \phi_{max})$ where $\phi_{max}$ is `priorParameters[1]`.
- For scale $\lambda$: $\text{Inv-Gamma}(a_\lambda, b_\lambda)$ where $a_\lambda$=`priorParameters[2]`, $b_\lambda$=`priorParameters[3]`.

Sampling procedure:
1. Draw shape $k$:
   $$k \sim U(0, \phi_{max})$$
2. Draw scale $\lambda$:
   $$\lambda \sim \text{Inv-Gamma}(a_\lambda, b_\lambda)$$
   (Equivalent to $1/\lambda \sim \text{Gamma}(a_\lambda, b_\lambda)$).

The pair $(k, \lambda)$ is one draw from $G_0$.
