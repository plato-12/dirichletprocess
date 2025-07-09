# Multivariate Normal Likelihood Function $F(\mathbf{y} | \boldsymbol{\mu}, \boldsymbol{\Sigma})$

## Mathematical Formula

For a $d$-dimensional data point $\mathbf{y}$, given a mean vector $\boldsymbol{\mu}$ ($d \times 1$) and a covariance matrix $\boldsymbol{\Sigma}$ ($d \times d$, symmetric positive definite), the likelihood is:

$$F(\mathbf{y} | \boldsymbol{\mu}, \boldsymbol{\Sigma}) = (2\pi)^{-d/2} |\boldsymbol{\Sigma}|^{-1/2} \exp\left(-\frac{1}{2}(\mathbf{y} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{y} - \boldsymbol{\mu})\right)$$

Where $|\cdot|$ denotes the determinant and $T$ denotes transpose.

## Implementation Considerations

1. Log-likelihood is numerically preferred:
   $$\log F(\mathbf{y} | \boldsymbol{\mu}, \boldsymbol{\Sigma}) = -\frac{d}{2}\log(2\pi) - \frac{1}{2}\log|\boldsymbol{\Sigma}| - \frac{1}{2}(\mathbf{y} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{y} - \boldsymbol{\mu})$$

2. The term $(\mathbf{y} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{y} - \boldsymbol{\mu})$ is the Mahalanobis distance squared.
3. Efficient computation of $\log|\boldsymbol{\Sigma}|$ and $\boldsymbol{\Sigma}^{-1}$ (or solving $\boldsymbol{\Sigma}\mathbf{x} = (\mathbf{y} - \boldsymbol{\mu})$ for $\mathbf{x}$ and then computing $\mathbf{x}^T(\mathbf{y} - \boldsymbol{\mu})$) is crucial. Often, the Cholesky decomposition of $\boldsymbol{\Sigma} = LL^T$ is used:
   - $\log|\boldsymbol{\Sigma}| = 2 \sum_j \log(L_{jj})$
   - The quadratic form can be computed by solving $L\mathbf{z} = (\mathbf{y} - \boldsymbol{\mu})$ for $\mathbf{z}$, then the quadratic form is $\mathbf{z}^T\mathbf{z}$.
4. Ensure $\boldsymbol{\Sigma}$ is symmetric and positive definite.
