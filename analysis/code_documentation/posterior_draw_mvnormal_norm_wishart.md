# Posterior Draw: MVN Parameters (Normal-Wishart Conjugacy)

Given data $D_k = \{\mathbf{y}_j : c_j=k\}$ for cluster $k$, $n_k = |D_k|$, sample mean $\bar{\mathbf{y}}_k = \frac{1}{n_k}\sum \mathbf{y}_j$, scatter matrix $\mathbf{S}_{scat,k} = \sum (\mathbf{y}_j - \bar{\mathbf{y}}_k)(\mathbf{y}_j - \bar{\mathbf{y}}_k)^T$.

Prior hyperparameters ($G_0$): $\boldsymbol{m}_0, \kappa_0, \nu_0, \mathbf{S}_0$ (scale matrix for precision $\boldsymbol{\Lambda}$).

Posterior hyperparameters for $(\boldsymbol{\mu}_k, \boldsymbol{\Lambda}_k = \boldsymbol{\Sigma}_k^{-1})$:
$$\kappa_n = \kappa_0 + n_k$$
$$\boldsymbol{m}_n = \frac{\kappa_0 \boldsymbol{m}_0 + n_k \bar{\mathbf{y}}_k}{\kappa_n}$$
$$\nu_n = \nu_0 + n_k$$
$$\mathbf{S}_n = \mathbf{S}_0 + \mathbf{S}_{scat,k} + \frac{\kappa_0 n_k}{\kappa_0 + n_k}(\bar{\mathbf{y}}_k - \boldsymbol{m}_0)(\bar{\mathbf{y}}_k - \boldsymbol{m}_0)^T$$

Sampling procedure:
1. Draw precision $\boldsymbol{\Lambda}_k | D_k, G_0 \sim \text{Wishart}(\nu_n, \mathbf{S}_n)$.
2. Draw mean $\boldsymbol{\mu}_k | \boldsymbol{\Lambda}_k, D_k, G_0 \sim N_d(\boldsymbol{m}_n, (\kappa_n \boldsymbol{\Lambda}_k)^{-1})$.
