# Posterior Draw: Weibull Parameters (Non-Conjugate via Metropolis-Hastings with Gibbs Step)

For a Weibull likelihood $F(y | k, \lambda)$ (shape $k$, scale $\lambda$) with potentially non-conjugate priors (e.g., $U(0, \phi_{max})$ for $k$, $\text{Inv-Gamma}(a_\lambda,b_\lambda)$ for $\lambda$).

Given data $D_k = \{y_j : c_j=k\}$ for cluster $k$.
Target posterior density:
$$p(k, \lambda | D_k, G_0) \propto \left( \prod_{y_j \in D_k} F(y_j | k, \lambda) \right) p_0(k) p_0(\lambda)$$

The `dirichletprocess` package often uses a Metropolis-within-Gibbs approach or a specialized MH for Weibull:
1. Sample shape $k'$ from a proposal $q(k'|k^{(t)})$ (e.g., random walk on $k$ or $\log k$).
2. Given $k'$, sample scale $\lambda'$ from its full conditional posterior $p(\lambda | k', D_k, G_0)$. If the prior $p_0(\lambda)$ is Inv-Gamma (i.e., $1/\lambda \sim \text{Gamma}(a_\lambda, b_\lambda)$), then the full conditional for $1/\lambda$ is:
   $$1/\lambda | k', D_k, G_0 \sim \text{Gamma}\left(a_\lambda + n_k, b_\lambda + \sum_{y_j \in D_k} (y_j)^{k'}\right)$$
   (The R code for `PosteriorDraw.weibull` samples $1/\lambda$ from $\text{Gamma}(n_k + a_0, \sum y_j^{k'} + b_0)$ implying $a_0, b_0$ are prior params for $1/\lambda$. The $(y_j/M_{scale})^{k'}$ form is more general if a scaling factor $M_{scale}$ were involved, but typically not for Weibull scale parameter $\lambda$ itself in this context).
3. Calculate acceptance ratio for $k'$ based on the marginal likelihood of $k'$ (integrating out $\lambda$) or using the joint proposal $(k', \lambda')$. The R code `MetropolisHastings.weibull` uses the joint proposal where $\lambda'$ is drawn from its conditional given $k'$.
