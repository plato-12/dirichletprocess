# Weibull Likelihood Function $F(y | k, \lambda)$

## Mathematical Formula

The Weibull distribution is often parameterized by a shape parameter $k > 0$ and a scale parameter $\lambda > 0$. The probability density function is:

$$F(y | k, \lambda) = \frac{k}{\lambda} \left(\frac{y}{\lambda}\right)^{k-1} \exp\left(-\left(\frac{y}{\lambda}\right)^k\right)$$

for $y \ge 0$.

The R `dweibull(x, shape, scale)` function uses this parameterization directly.

## Implementation Considerations

1. Ensure $y \ge 0$, $k > 0$, and $\lambda > 0$.
2. Log-likelihood:
   $$\log F(y | k, \lambda) = \log(k) - \log(\lambda) + (k-1)(\log(y) - \log(\lambda)) - \left(\frac{y}{\lambda}\right)^k$$
3. Be careful with $y=0$: If $k=1$ (Exponential), $f(0|1,\lambda)=1/\lambda$. If $k>1$, $f(0|k,\lambda)=0$. If $0<k<1$, $f(0|k,\lambda)=\infty$.
