# Posterior Draw: Gaussian Parameters (Normal-Inverse Gamma Conjugacy)

Given data $D_k = \{y_j : c_j=k\}$ assigned to cluster $k$, with $n_k = |D_k|$, sample mean $\bar{y}_k = \frac{1}{n_k}\sum_{y_j \in D_k} y_j$, and sum of squared deviations $S_k = \sum_{y_j \in D_k} (y_j - \bar{y}_k)^2$.

Prior hyperparameters (from $G_0$): $\mu_0, \kappa_0, a_0, b_0$.

Posterior hyperparameters:
$$\kappa_n = \kappa_0 + n_k$$
$$\mu_n = \frac{\kappa_0 \mu_0 + n_k \bar{y}_k}{\kappa_n}$$
$$a_n = a_0 + \frac{n_k}{2}$$
$$b_n = b_0 + \frac{1}{2}S_k + \frac{\kappa_0 n_k (\bar{y}_k - \mu_0)^2}{2(\kappa_0 + n_k)}$$

Sampling procedure for cluster parameters $(\mu_k, \sigma_k^2)$:
1. Draw variance $\sigma_k^2$ from the posterior Inverse Gamma distribution:
   $$\sigma_k^2 | D_k, G_0 \sim \text{Inv-Gamma}(a_n, b_n)$$
2. Draw mean $\mu_k$ conditional on $\sigma_k^2$ from the posterior Normal distribution:
   $$\mu_k | \sigma_k^2, D_k, G_0 \sim N(\mu_n, \sigma_k^2 / \kappa_n)$$
