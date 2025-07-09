#include "../inst/include/weibull_mixing.h"
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

namespace dirichletprocess {

// Constructor
WeibullMixing::WeibullMixing(double phi, double alpha0, double beta0,
                             double hyper_a1, double hyper_a2,
                             double hyper_b1, double hyper_b2,
                             double mh_step_alpha, int mh_draws)
  : phi(phi), alpha0(alpha0), beta0(beta0),
    hyper_a1(hyper_a1), hyper_a2(hyper_a2),
    hyper_b1(hyper_b1), hyper_b2(hyper_b2),
    mh_step_alpha(mh_step_alpha), mh_draws(mh_draws) {}

// Optimized log likelihood implementation
double WeibullMixing::log_likelihood(const arma::vec& data_point,
                                     const arma::vec& params) const {
  double x = data_point[0];
  double alpha = params[0];
  double lambda = params[1];

  // Check bounds
  if (x <= 0 || alpha <= 0 || lambda <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Optimized Weibull log-likelihood
  // Use log(x) to avoid expensive pow() call
  double log_x = std::log(x);
  double log_lik = -std::log(lambda) + std::log(alpha) +
    (alpha - 1.0) * log_x -
    std::exp(alpha * log_x) / lambda;

  return log_lik;
}

// Prior draw
arma::vec WeibullMixing::prior_draw() const {
  arma::vec params(2);

  // alpha ~ Uniform(0, phi)
  params[0] = R::runif(0, phi);

  // lambda = 1/Gamma(alpha0, beta0)
  double gamma_draw = R::rgamma(alpha0, 1.0 / beta0);
  params[1] = 1.0 / std::max(1e-10, gamma_draw);

  return params;
}

// Optimized posterior draw using Metropolis-Hastings
arma::vec WeibullMixing::posterior_draw(const arma::mat& cluster_data,
                                        const arma::vec& prior_params) const {
  int n = cluster_data.n_rows;

  // Handle empty cluster
  if (n == 0) {
    return prior_draw();
  }

  // Pre-compute log(x) for all data points to avoid repeated calculations
  arma::vec log_x(n);
  for (int i = 0; i < n; ++i) {
    double xi = cluster_data(i, 0);
    log_x[i] = (xi > 0) ? std::log(xi) : -std::numeric_limits<double>::infinity();
  }

  // Initialize with prior draw
  arma::vec current_params = prior_draw();

  // Run Metropolis-Hastings with Gibbs update for lambda
  for (int iter = 0; iter < mh_draws; ++iter) {
    double alpha_current = current_params[0];

    // Propose new alpha
    double alpha_prop = std::abs(alpha_current + mh_step_alpha * R::rnorm(0, 1.7));
    if (alpha_prop > phi) {
      alpha_prop = phi;
    }

    // Efficiently compute sum(x^alpha) using pre-computed log(x)
    double sum_x_alpha_current = 0.0;
    double sum_x_alpha_prop = 0.0;

    for (int i = 0; i < n; ++i) {
      if (std::isfinite(log_x[i])) {
        sum_x_alpha_current += std::exp(alpha_current * log_x[i]);
        sum_x_alpha_prop += std::exp(alpha_prop * log_x[i]);
      }
    }

    // Sample lambda values using conjugate posteriors
    double shape_post = n + alpha0;
    double rate_post_current = sum_x_alpha_current + beta0;
    double rate_post_prop = sum_x_alpha_prop + beta0;

    double gamma_current = R::rgamma(shape_post, 1.0 / rate_post_current);
    double lambda_current = 1.0 / std::max(1e-10, gamma_current);

    double gamma_prop = R::rgamma(shape_post, 1.0 / rate_post_prop);
    double lambda_prop = 1.0 / std::max(1e-10, gamma_prop);

    // Compute log likelihoods efficiently
    double log_lik_current = 0.0;
    double log_lik_prop = 0.0;

    // Pre-compute constants
    double log_lambda_current = std::log(lambda_current);
    double log_lambda_prop = std::log(lambda_prop);
    double log_alpha_current = std::log(alpha_current);
    double log_alpha_prop = std::log(alpha_prop);

    for (int i = 0; i < n; ++i) {
      if (std::isfinite(log_x[i])) {
        // Current parameters
        log_lik_current += log_alpha_current - log_lambda_current +
          (alpha_current - 1.0) * log_x[i] -
          std::exp(alpha_current * log_x[i]) / lambda_current;

        // Proposed parameters
        log_lik_prop += log_alpha_prop - log_lambda_prop +
          (alpha_prop - 1.0) * log_x[i] -
          std::exp(alpha_prop * log_x[i]) / lambda_prop;
      }
    }

    // Compute log priors for alpha only
    double log_prior_alpha_current = (alpha_current > 0 && alpha_current <= phi) ?
    -std::log(phi) : -std::numeric_limits<double>::infinity();
    double log_prior_alpha_prop = (alpha_prop > 0 && alpha_prop <= phi) ?
    -std::log(phi) : -std::numeric_limits<double>::infinity();

    // Accept/reject
    double log_ratio = (log_lik_prop + log_prior_alpha_prop) -
    (log_lik_current + log_prior_alpha_current);

    if (!std::isfinite(log_ratio)) {
      log_ratio = -std::numeric_limits<double>::infinity();
    }

    double accept_prob = std::min(1.0, std::exp(log_ratio));

    if (R::runif(0, 1) < accept_prob) {
      current_params[0] = alpha_prop;
      current_params[1] = lambda_prop;
    } else {
      current_params[0] = alpha_current;
      current_params[1] = lambda_current;
    }
  }

  return current_params;
}

// Log prior density
double WeibullMixing::log_prior_density(const arma::vec& params) const {
  double alpha = params[0];
  double lambda = params[1];

  // Check bounds
  if (alpha <= 0 || alpha > phi || lambda <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Log prior for alpha: log(1/phi) = -log(phi)
  double log_prior_alpha = -std::log(phi);

  // Log prior for lambda: Inverse-Gamma(alpha0, beta0)
  double log_prior_lambda = alpha0 * std::log(beta0) - std::lgamma(alpha0) -
    (alpha0 + 1.0) * std::log(lambda) - beta0 / lambda;

  return log_prior_alpha + log_prior_lambda;
}

// MH parameter proposal
arma::vec WeibullMixing::mh_parameter_proposal(const arma::vec& current_params) const {
  arma::vec proposed = current_params;

  // Propose new alpha
  double alpha_current = current_params[0];
  double alpha_proposal = std::abs(alpha_current + mh_step_alpha * R::rnorm(0, 1.7));

  // Ensure alpha stays within bounds
  if (alpha_proposal > phi) {
    alpha_proposal = phi;
  }

  proposed[0] = alpha_proposal;
  // Lambda will be updated analytically

  return proposed;
}

// Update hyperparameters
void WeibullMixing::update_hyperparameters(const std::vector<arma::vec>& all_params) {
  int K = all_params.size();
  if (K == 0) return;

  // Find maximum alpha and sum of inverse lambdas
  double max_alpha = 0.0;
  double sum_inv_lambda = 0.0;

  for (const auto& params : all_params) {
    max_alpha = std::max(max_alpha, params[0]);
    if (params[1] > 1e-10) {
      sum_inv_lambda += 1.0 / params[1];
    }
  }

  // Update phi using Pareto posterior
  double xm = std::max(max_alpha, hyper_a1);
  double shape = hyper_a2 + K;
  double U = R::runif(0, 1);
  phi = qpareto(U, xm, shape);

  // Update beta0 using Gamma posterior
  double post_shape = hyper_b1 + 2 * K;
  double post_rate = hyper_b2 + sum_inv_lambda;
  beta0 = R::rgamma(post_shape, 1.0 / post_rate);
}

// Pareto quantile function
double WeibullMixing::qpareto(double p, double xm, double alpha) const {
  if (p <= 0 || p >= 1) {
    Rcpp::stop("p must be in (0,1) for qpareto");
  }
  return xm * std::pow(1 - p, -1.0 / alpha);
}

} // namespace dirichletprocess
