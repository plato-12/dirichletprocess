# Gaussian Likelihood Function $F(y | \mu, \sigma^2)$

## Mathematical Formula

The likelihood for a data point $y$ given mean $\mu$ and variance $\sigma^2$ (or standard deviation $\sigma$) is:

$$F(y | \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(y - \mu)^2}{2\sigma^2}\right) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(y - \mu)^2}{2\sigma^2}\right)$$

## Implementation Considerations

1. For numerical stability, compute log-likelihood first:
   $$\log F(y | \mu, \sigma^2) = -\frac{1}{2}\log(2\pi\sigma^2) - \frac{(y - \mu)^2}{2\sigma^2} = -\log(\sigma) - \frac{1}{2}\log(2\pi) - \frac{(y - \mu)^2}{2\sigma^2}$$

2. Vectorization: Compute for multiple data points $y_i$ simultaneously given a single $(\mu, \sigma^2)$ or for a single $y_i$ against multiple $(\mu_k, \sigma_k^2)$.

3. Ensure $\sigma^2 > 0$ (or $\sigma > 0$). Handle cases where $\sigma$ is very small to avoid overflow/underflow.
