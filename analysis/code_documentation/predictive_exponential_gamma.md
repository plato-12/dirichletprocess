# Predictive Likelihood: New Point from Exponential with Gamma Prior

For an Exponential likelihood $F(y | \lambda_{rate})$ and a Gamma prior $G_0 = \text{Gamma}(\lambda_{rate} | a_0, b_0)$ on the rate.

The predictive likelihood for a new data point $y_{new}$ is:
$$p(y_{new} | G_0) = \int_0^\infty (\lambda_{rate} e^{-\lambda_{rate} y_{new}}) \frac{b_0^{a_0}}{\Gamma(a_0)} \lambda_{rate}^{a_0-1} e^{-b_0 \lambda_{rate}} d\lambda_{rate}$$
$$= \frac{b_0^{a_0}}{\Gamma(a_0)} \int_0^\infty \lambda_{rate}^{a_0} e^{-(b_0+y_{new})\lambda_{rate}} d\lambda_{rate}$$
This is related to the Gamma function. The integral is $\frac{\Gamma(a_0+1)}{(b_0+y_{new})^{a_0+1}}$.
So, $$p(y_{new} | G_0) = \frac{b_0^{a_0}}{\Gamma(a_0)} \frac{\Gamma(a_0+1)}{(b_0+y_{new})^{a_0+1}} = \frac{a_0 b_0^{a_0}}{(b_0+y_{new})^{a_0+1}}$$
This is the PDF of a Lomax (Pareto Type II) distribution. The R code `Predictive.exponential` calculates $a_n = a_0+1, b_n = b_0+y_{new}$ and then uses $\frac{\Gamma(a_n)}{\Gamma(a_0)} \frac{b_0^{a_0}}{b_n^{a_n}}$.
