# Beta Likelihood Function $F(y | \mu, \tau, M_{max})$

## Mathematical Formula

The `dirichletprocess` package often parameterizes the Beta distribution for data $y \in [0, M_{max}]$ using a mean $\mu_y \in [0, M_{max}]$ and a precision parameter $\tau > 0$. Let $x = y/M_{max}$ be the scaled data point in $[0,1]$, with mean $\mu_x = \mu_y/M_{max}$. The standard Beta distribution for $x$ is $\text{Beta}(x | a, b)$. The parameters $a$ and $b$ are related to $\mu_x$ and $\tau$ by:

$$a = \mu_x \tau$$
$$b = (1 - \mu_x) \tau$$

The likelihood for $y$ is then given by the transformed density:

$$F(y | \mu_y, \tau, M_{max}) = \frac{1}{M_{max}} \text{Beta}\left(\frac{y}{M_{max}} \Big| a = \frac{\mu_y}{M_{max}}\tau, b = \left(1 - \frac{\mu_y}{M_{max}}\right)\tau \right)$$
$$= \frac{1}{M_{max}} \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \left(\frac{y}{M_{max}}\right)^{a-1} \left(1 - \frac{y}{M_{max}}\right)^{b-1}$$

Where $0 \le y \le M_{max}$, $0 < \mu_y < M_{max}$, and $\tau > 0$.

## Implementation Considerations

1. Ensure $a > 0$ and $b > 0$. This implies $\mu_y \in (0, M_{max})$ and $\tau > 0$.
2. Handle edge cases: If $y=0$ or $y=M_{max}$, the likelihood can be 0 or infinity if $a<1$ or $b<1$. Standard library functions for Beta PDF usually handle this.
3. Log-likelihood is preferred for products: $\log F(y | \dots) = -\log(M_{max}) + \log\Gamma(a+b) - \log\Gamma(a) - \log\Gamma(b) + (a-1)\log(y/M_{max}) + (b-1)\log(1-y/M_{max})$.
