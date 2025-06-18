#include "../inst/include/beta_mixing.h"
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

namespace dirichletprocess {

// Constructor
BetaMixing::BetaMixing(double alpha0, double beta0, double maxT)
  : alpha0(alpha0), beta0(beta0), maxT(maxT) {}

// Log likelihood implementation
double BetaMixing::log_likelihood(const arma::vec& data_point,
                                  const arma::vec& params) const {
  double mu = params[0];
  double tau = params[1];

  // Convert to standard Beta parameters
  double a = (mu * tau) / maxT;
  double b = (1.0 - mu/maxT) * tau;

  double x = data_point[0];

  // Check bounds
  if (x < 0 || x > maxT || a <= 0 || b <= 0 || tau <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Log likelihood: log(1/maxT) + log(dbeta(x/maxT, a, b))
  return std::log(1.0/maxT) + R::dbeta(x/maxT, a, b, 1);
}

// Posterior draw - for non-conjugate, we typically use MH or return prior draw
arma::vec BetaMixing::posterior_draw(const arma::mat& cluster_data,
                                     const arma::vec& prior_params) const {
  // For non-conjugate Beta, we would typically use Metropolis-Hastings
  // For now, return a simple approximation based on method of moments
  if (cluster_data.n_rows == 0) {
    return prior_draw();
  }

  arma::vec x = cluster_data.col(0);

  // Method of moments estimates
  double x_mean = arma::mean(x) / maxT;  // Normalize to [0,1]
  double x_var = arma::var(x) / (maxT * maxT);

  // Prevent division by zero
  if (x_var < 1e-10) {
    x_var = 1e-10;
  }

  // Calculate tau (precision) from variance
  double tau_est = x_mean * (1 - x_mean) / x_var - 1;
  tau_est = std::max(0.1, tau_est);  // Ensure positive

  // mu is mean * maxT
  double mu_est = x_mean * maxT;

  // Add some noise for MCMC
  arma::vec params(2);
  params[0] = std::max(0.01, std::min(maxT - 0.01,
                                      mu_est + 0.1 * R::rnorm(0, 1)));
  params[1] = std::max(0.1, tau_est + 0.1 * R::rnorm(0, 1));

  return params;
}

// Prior draw
arma::vec BetaMixing::prior_draw() const {
  arma::vec params(2);

  // mu ~ Uniform(0, maxT)
  params[0] = R::runif(0, maxT);

  // tau ~ InverseGamma(alpha0, beta0), so 1/tau ~ Gamma(alpha0, beta0)
  double gamma_draw = R::rgamma(alpha0, 1.0/beta0);  // R uses scale, not rate
  params[1] = 1.0 / std::max(1e-10, gamma_draw);     // tau = 1/gamma

  return params;
}

} // namespace dirichletprocess
