# Prior Draw: Gaussian Parameters from Normal-Inverse Gamma ($G_0$)

For a Gaussian likelihood $N(y | \mu, \sigma^2)$, a common conjugate base measure $G_0$ is the Normal-Inverse Gamma distribution for $(\mu, \sigma^2)$.

Hyperparameters for $G_0$: $\mu_0, \kappa_0 > 0, a_0 > 0, b_0 > 0$.

Sampling procedure:
1. Draw variance $\sigma^2$ from an Inverse Gamma distribution:
   $$\sigma^2 \sim \text{Inv-Gamma}(a_0, b_0)$$
   (This is equivalent to drawing precision $\lambda = 1/\sigma^2 \sim \text{Gamma}(a_0, b_0)$).
2. Draw mean $\mu$ conditional on the drawn $\sigma^2$ from a Normal distribution:
   $$\mu | \sigma^2 \sim N(\mu_0, \sigma^2 / \kappa_0)$$

The pair $(\mu, \sigma^2)$ is one draw from $G_0$.
