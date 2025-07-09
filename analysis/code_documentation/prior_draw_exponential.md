# Prior Draw: Exponential Parameter (rate $\lambda_{rate}$) from $G_0$

For an Exponential likelihood $F(y | \lambda_{rate})$, a common conjugate base measure $G_0$ for the rate parameter $\lambda_{rate}$ is a Gamma distribution.

Hyperparameters for $G_0$: shape $a_0 > 0$, rate $b_0 > 0$.

Sampling procedure:
1. Draw rate $\lambda_{rate}$ from a Gamma distribution:
   $$\lambda_{rate} \sim \text{Gamma}(a_0, b_0)$$

This $\lambda_{rate}$ is one draw from $G_0$.
