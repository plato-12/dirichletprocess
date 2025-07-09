// src/gaussian_mixing.cpp
#include <RcppArmadillo.h>
#include "../inst/include/gaussian_mixing.h"

namespace dirichletprocess {

GaussianMixing::GaussianMixing(double mu0, double kappa0,
                               double alpha0, double beta0)
  : mu0(mu0), kappa0(kappa0), alpha0(alpha0), beta0(beta0) {}

double GaussianMixing::log_likelihood(const arma::vec& data_point,
                                      const arma::vec& params) const {
  if (params.n_elem < 2) {
    return -std::numeric_limits<double>::infinity();
  }

  double mean = params[0];
  double variance = params[1];

  // Ensure variance is positive
  if (variance <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Compute log likelihood with numerical stability
  double standardized = (data_point[0] - mean) / std::sqrt(variance);
  double log_lik = -0.5 * (std::log(2 * M_PI) + std::log(variance) + standardized * standardized);

  // Check for numerical issues
  if (!std::isfinite(log_lik)) {
    return -std::numeric_limits<double>::infinity();
  }

  return log_lik;
}

arma::vec GaussianMixing::posterior_draw(const arma::mat& cluster_data,
                                         const arma::vec& prior_params) const {
  int n = cluster_data.n_rows;
  double data_mean = arma::mean(cluster_data.col(0));
  double data_var = arma::var(cluster_data.col(0));

  // Posterior parameters (Normal-Inverse-Gamma)
  double kappa_n = kappa0 + n;
  double mu_n = (kappa0 * mu0 + n * data_mean) / kappa_n;
  double alpha_n = alpha0 + n / 2.0;
  double beta_n = beta0 + 0.5 * n * data_var +
    0.5 * kappa0 * n * std::pow(data_mean - mu0, 2) / kappa_n;

  // Sample variance from Inverse-Gamma
  double variance = 1.0 / R::rgamma(alpha_n, 1.0 / beta_n);

  // Sample mean from Normal
  double mean = R::rnorm(mu_n, std::sqrt(variance / kappa_n));

  arma::vec params(2);
  params[0] = mean;
  params[1] = variance;

  return params;
}

arma::vec GaussianMixing::prior_draw() const {
  double variance = 1.0 / R::rgamma(alpha0, 1.0 / beta0);
  double mean = R::rnorm(mu0, std::sqrt(variance / kappa0));

  arma::vec params(2);
  params[0] = mean;
  params[1] = variance;

  return params;
}

} // namespace dirichletprocess
