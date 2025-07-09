# Exponential Likelihood Function $F(y | \lambda_{rate})$

## Mathematical Formula

The Exponential distribution is parameterized by a rate parameter $\lambda_{rate} > 0$. The probability density function is:

$$F(y | \lambda_{rate}) = \lambda_{rate} e^{-\lambda_{rate} y}$$

for $y \ge 0$.

## Implementation Considerations

1. Ensure $y \ge 0$ and $\lambda_{rate} > 0$.
2. Log-likelihood: $\log F(y | \lambda_{rate}) = \log(\lambda_{rate}) - \lambda_{rate} y$.
