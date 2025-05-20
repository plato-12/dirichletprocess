# Posterior Draw: Exponential Parameter (Gamma Conjugacy)

For an Exponential likelihood $F(y | \lambda_{rate})$ and a Gamma prior $G_0 = \text{Gamma}(\lambda_{rate} | a_0, b_0)$ on the rate parameter.

Given data $D_k = \{y_j : c_j=k\}$ for cluster $k$, with $n_k = |D_k|$ and sum $S_k = \sum_{y_j \in D_k} y_j$.

Posterior hyperparameters for $\lambda_{rate,k}$:
$$a_n = a_0 + n_k$$
$$b_n = b_0 + S_k$$

Sampling procedure for cluster rate parameter $\lambda_{rate,k}$:
1. Draw $\lambda_{rate,k}$ from the posterior Gamma distribution:
   $$\lambda_{rate,k} | D_k, G_0 \sim \text{Gamma}(a_n, b_n)$$
