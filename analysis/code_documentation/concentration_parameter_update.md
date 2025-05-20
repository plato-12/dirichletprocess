# Concentration Parameter $\alpha$ Update (Escobar & West, 1995)

The concentration parameter $\alpha$ controls the expected number of clusters. It is typically assigned a Gamma prior: $\alpha \sim \text{Gamma}(a_\alpha, b_\alpha)$.

Given $K$ current clusters among $N$ data points, the update uses an auxiliary variable $\eta$. Following Escobar and West (1995), and as implemented in `dirichletprocess` (see `update_concentration` R function):

1. Draw $\eta | \alpha, N \sim \text{Beta}(\alpha + 1, N)$.
2. The posterior for $\alpha$ is then a mixture of two Gamma distributions:
   $$p(\alpha | K, N, \eta, a_\alpha, b_\alpha) = \pi_\eta \cdot \text{Gamma}(\alpha | a_\alpha+K, b_\alpha-\log(\eta)) + (1-\pi_\eta) \cdot \text{Gamma}(\alpha | a_\alpha+K-1, b_\alpha-\log(\eta))$$
   where the mixture probability $\pi_\eta$ is given by:
   $$\pi_\eta = \frac{a_\alpha+K-1}{(a_\alpha+K-1) + N(b_\alpha-\log(\eta))}$$
   (This assumes $a_\alpha+K-1 > 0$. If $K=0$, the second Gamma component might use shape $a_\alpha$. The R code handles cases where $a_\alpha+K-1 \le 0$ by adjusting which component is chosen or falling back to the prior).

A draw is made from this mixture distribution to get the new $\alpha$.
