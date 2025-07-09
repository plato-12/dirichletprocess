# Predictive Likelihood: New Point from Gaussian with Normal-Inverse Gamma Prior

The predictive likelihood for a new data point $y_{new}$, given $G_0 \sim \text{Normal-Inverse Gamma}(\mu_0, \kappa_0, a_0, b_0)$, is the density of $y_{new}$ under a non-standardized Student's t-distribution:

$$p(y_{new} | G_0) = \int N(y_{new} | \mu, \sigma^2) p(\mu, \sigma^2 | G_0) d\mu d\sigma^2$$
This evaluates to a density proportional to:
$$p(y_{new} | G_0) \propto \left(1 + \frac{\kappa_0(y_{new}-\mu_0)^2}{(\kappa_0+1)2b_0}\right)^{-(a_0+1/2)}$$
The exact form used in `Predictive.normal` involves calculating posterior parameters $(a_n, b_n, \kappa_n)$ as if $y_{new}$ was observed, and then using a ratio of Gamma functions and other terms:
$$p(y_{new}|G_0) = \frac{\Gamma(a_n)}{\Gamma(a_0)} \frac{b_0^{a_0}}{b_n^{a_n}} \sqrt{\frac{\kappa_0}{\kappa_n (2\pi)}}$$ (The $2\pi$ factor might be handled differently or cancel in practice within the CRP probabilities if not careful with exact density vs. proportional to).
The `dirichletprocess` R code for `Predictive.normal` is: $\frac{\Gamma(a_n)}{\Gamma(a_0)} \frac{b_0^{a_0}}{b_n^{a_n}} \sqrt{\frac{\kappa_0}{\kappa_n}}$, where $a_n, b_n, \kappa_n$ are posterior parameters after observing $y_{new}$. This form is common for the marginal likelihood used in Bayes Factors or model comparison, and is proportional to the Student-t PDF.
