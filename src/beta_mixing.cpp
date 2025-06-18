#include "../inst/include/beta_mixing.h"
#include <RcppArmadillo.h>
#include <cmath>

namespace dirichletprocess {

double BetaMixing::log_likelihood(const arma::vec& data_point,
                                  const arma::vec& params) const {
  double x = data_point[0];
  double mu = params[0];
  double tau = params[1];

  // Transform to standard Beta parameters
  double a = (mu * tau) / maxT;
  double b = (1.0 - mu/maxT) * tau;

  // Bounds checking
  if (x < 0.0 || x > maxT) return -std::numeric_limits<double>::infinity();
  if (a <= 0.0 || b <= 0.0) return -std::numeric_limits<double>::infinity();
  if (mu < 0.0 || mu > maxT) return -std::numeric_limits<double>::infinity();
  if (tau <= 0.0) return -std::numeric_limits<double>::infinity();

  // Scaled Beta log-likelihood
  double x_scaled = x / maxT;
  double log_lik = (a - 1.0) * std::log(x_scaled) +
    (b - 1.0) * std::log(1.0 - x_scaled) +
    lgamma(a + b) - lgamma(a) - lgamma(b) -
    std::log(maxT);

  return log_lik;
}

arma::vec BetaMixing::posterior_draw(const arma::mat& cluster_data,
                                     const arma::vec& prior_params) const {
  int n = cluster_data.n_rows;

  // Handle empty cluster case
  if (n == 0) {
    return prior_draw();
  }

  // Since Beta has no conjugate prior, use Metropolis-Hastings
  // Start from method of moments estimate
  arma::vec start_params = method_of_moments_estimate(cluster_data.col(0));

  // Default step sizes
  arma::vec step_sizes = {1.0, 1.0};

  // Run Metropolis-Hastings
  return metropolis_hastings_step(cluster_data, start_params, step_sizes, 250);
}

arma::vec BetaMixing::prior_draw() const {
  arma::vec params(2);

  // Draw mu ~ Uniform(0, maxT)
  params[0] = R::runif(0.0, maxT);

  // Draw nu ~ InverseGamma(alpha0, beta0)
  // Since nu = 1/tau and tau ~ Gamma(alpha0, beta0)
  double tau = R::rgamma(alpha0, 1.0/beta0);
  params[1] = tau;

  return params;
}

arma::vec BetaMixing::method_of_moments_estimate(const arma::vec& x) const {
  double x_mean = arma::mean(x);
  double x_var = arma::var(x);

  // Prevent division by zero
  if (x_var < 1e-10) {
    x_var = 1e-10;
  }

  // Method of moments for Beta distribution
  double x_mean_scaled = x_mean / maxT;
  double x_var_scaled = x_var / (maxT * maxT);

  // Estimate parameters
  double common = x_mean_scaled * (1.0 - x_mean_scaled) / x_var_scaled - 1.0;
  common = std::max(0.1, common);  // Ensure positive

  arma::vec params(2);
  params[0] = x_mean;  // mu
  params[1] = common;  // tau

  return params;
}

arma::vec BetaMixing::metropolis_hastings_step(const arma::mat& cluster_data,
                                               const arma::vec& current_params,
                                               const arma::vec& step_sizes,
                                               int n_draws) const {
  int n = cluster_data.n_rows;
  arma::vec x = cluster_data.col(0);

  arma::vec params = current_params;
  int accepted = 0;

  for (int i = 0; i < n_draws; ++i) {
    // Propose new parameters
    arma::vec proposal = params;
    proposal[0] += step_sizes[0] * R::rnorm(0.0, 2.4);
    proposal[1] = std::abs(params[1] + step_sizes[1] * R::rnorm(0.0, 2.4));

    // Check bounds for mu
    if (proposal[0] < 0.0 || proposal[0] > maxT) {
      continue;
    }

    // Calculate log likelihood ratio
    double log_ratio = 0.0;
    for (int j = 0; j < n; ++j) {
      arma::vec data_point = x.row(j).t();
      log_ratio += log_likelihood(data_point, proposal) -
        log_likelihood(data_point, params);
    }

    // Add prior ratio
    // Prior: mu ~ Uniform(0, maxT), tau ~ Gamma(alpha0, beta0)
    log_ratio += (alpha0 - 1.0) * (std::log(proposal[1]) - std::log(params[1]));
    log_ratio -= beta0 * (proposal[1] - params[1]);

    // Accept/reject
    if (std::log(R::runif(0.0, 1.0)) < log_ratio) {
      params = proposal;
      accepted++;
    }
  }

  return params;
}

} // namespace dirichletprocess
