# Prior Draw: Multivariate Normal Parameters from Normal-Wishart ($G_0$)

For an MVN likelihood $N_d(\mathbf{y} | \boldsymbol{\mu}, \boldsymbol{\Sigma})$, a conjugate base measure $G_0$ for $(\boldsymbol{\mu}, \boldsymbol{\Lambda} = \boldsymbol{\Sigma}^{-1})$ is Normal-Wishart.

Hyperparameters for $G_0$: mean vector $\boldsymbol{m}_0$, scalar $\kappa_0 > 0$, degrees of freedom $\nu_0 > d-1$, scale matrix $\mathbf{S}_0$ ($d \times d$, symmetric positive definite for $\boldsymbol{\Lambda}$).

Sampling procedure:
1. Draw precision matrix $\boldsymbol{\Lambda}$ from a Wishart distribution:
   $$\boldsymbol{\Lambda} \sim \text{Wishart}(\nu_0, \mathbf{S}_0)$$
   (Note: R's `rWishart(nu, V)` samples from Wishart with df $\nu$ and scale matrix $V$. If $\mathbf{S}_0$ is the scale matrix for $\boldsymbol{\Lambda}$, this is direct. Sometimes $\mathbf{S}_0$ is parameterized as $(E[\boldsymbol{\Lambda}])^{-1}/\nu_0$).

2. Draw mean vector $\boldsymbol{\mu}$ conditional on $\boldsymbol{\Lambda}$ from a Multivariate Normal distribution:
   $$\boldsymbol{\mu} | \boldsymbol{\Lambda} \sim N_d(\boldsymbol{m}_0, (\kappa_0 \boldsymbol{\Lambda})^{-1})$$

The pair $(\boldsymbol{\mu}, \boldsymbol{\Sigma} = \boldsymbol{\Lambda}^{-1})$ is one draw from $G_0$.
