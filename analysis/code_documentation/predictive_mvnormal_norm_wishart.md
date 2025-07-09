# Predictive Likelihood: New Point from MVN with Normal-Wishart Prior

For an MVN likelihood $N_d(\mathbf{y} | \boldsymbol{\mu}, \boldsymbol{\Sigma})$ and a Normal-Wishart prior $G_0$ on $(\boldsymbol{\mu}, \boldsymbol{\Lambda}=\boldsymbol{\Sigma}^{-1})$ with hyperparameters $(\boldsymbol{m}_0, \kappa_0, \nu_0, \mathbf{S}_0)$.

The predictive likelihood for a new data point $\mathbf{y}_{new}$ is the density of $\mathbf{y}_{new}$ under a multivariate Student's t-distribution:
$$p(\mathbf{y}_{new} | G_0) = \int N_d(\mathbf{y}_{new} | \boldsymbol{\mu}, \boldsymbol{\Lambda}^{-1}) p(\boldsymbol{\mu}, \boldsymbol{\Lambda} | G_0) d\boldsymbol{\mu} d\boldsymbol{\Lambda}$$
This evaluates to $t_{\nu_0-d+1}\left(\mathbf{y}_{new} \Big| \boldsymbol{m}_0, \frac{\mathbf{S}_0^{-1}(\kappa_0+1)}{\kappa_0(\nu_0-d+1)}\right)$ (scaled and shifted multivariate t).

The R code `Predictive.mvnormal` calculates posterior parameters $(\boldsymbol{m}_n, \kappa_n, \nu_n, \mathbf{S}_n)$ as if $\mathbf{y}_{new}$ were observed, and then uses a formula involving ratios of determinants and Gamma functions:
$$p(\mathbf{y}_{new} | G_0) = \pi^{-d/2} \frac{\Gamma_d(\nu_n/2)}{\Gamma_d(\nu_0/2)} \frac{|\mathbf{S}_0|^{\nu_0/2}}{|\mathbf{S}_n|^{\nu_n/2}} \left(\frac{\kappa_0}{\kappa_n}\right)^{d/2}$$
where $\Gamma_d(\cdot)$ is the multivariate Gamma function. This is proportional to the multivariate Student-t PDF.
