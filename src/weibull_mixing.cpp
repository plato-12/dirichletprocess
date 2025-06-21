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

// Log likelihood implementation
double WeibullMixing::log_likelihood(const arma::vec& data_point,
                                     const arma::vec& params) const {
  double x = data_point[0];
  double alpha = params[0];
  double lambda = params[1];

  // Check bounds
  if (x <= 0 || alpha <= 0 || lambda <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Weibull PDF: f(x|α,λ) = (α/λ) * (x/λ)^(α-1) * exp(-(x/λ)^α)
  // Log PDF: log(α) - log(λ) + (α-1)*[log(x) - log(λ)] - (x/λ)^α
  double log_lik = std::log(alpha) - std::log(lambda) +
    (alpha - 1.0) * (std::log(x) - std::log(lambda)) -
    std::pow(x / lambda, alpha);

  return log_lik;
}

// Prior draw
arma::vec WeibullMixing::prior_draw() const {
  arma::vec params(2);

  // alpha ~ Uniform(0, phi)
  params[0] = R::runif(0, phi);

  // 1/lambda ~ Gamma(alpha0, beta0), so lambda = 1/Gamma(alpha0, beta0)
  double gamma_draw = R::rgamma(alpha0, 1.0/beta0);  // R uses scale parameterization
  params[1] = 1.0 / std::max(1e-10, gamma_draw);

  return params;
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
  // p(lambda) = (beta0^alpha0 / Gamma(alpha0)) * lambda^(-alpha0-1) * exp(-beta0/lambda)
  double log_prior_lambda = alpha0 * std::log(beta0) - std::lgamma(alpha0) -
    (alpha0 + 1.0) * std::log(lambda) - beta0 / lambda;

  return log_prior_alpha + log_prior_lambda;
}

// MH parameter proposal
arma::vec WeibullMixing::mh_parameter_proposal(const arma::vec& current_params) const {
  arma::vec proposed = current_params;

  // Propose new alpha with reflecting boundary
  double alpha_current = current_params[0];
  double alpha_proposal = alpha_current + mh_step_alpha * R::rnorm(0, 1.7);

  // Reflect at boundaries
  if (alpha_proposal < 0) {
    alpha_proposal = -alpha_proposal;
  } else if (alpha_proposal > phi) {
    alpha_proposal = 2 * phi - alpha_proposal;
  }

  proposed[0] = alpha_proposal;
  // Lambda will be updated analytically given alpha

  return proposed;
}

// Posterior draw using Metropolis-Hastings
arma::vec WeibullMixing::posterior_draw(const arma::mat& cluster_data,
                                        const arma::vec& prior_params) const {
  int n = cluster_data.n_rows;

  // Handle empty cluster
  if (n == 0) {
    return prior_draw();
  }

  // Initialize with prior draw
  arma::vec current_params = prior_draw();
  double current_log_lik = 0.0;
  double current_log_prior = log_prior_density(current_params);

  // Calculate initial log likelihood
  for (int i = 0; i < n; ++i) {
    current_log_lik += log_likelihood(cluster_data.row(i).t(), current_params);
  }

  // Metropolis-Hastings iterations
  int accept_count = 0;

  for (int iter = 0; iter < mh_draws; ++iter) {
    // Propose new alpha
    arma::vec proposed_params = mh_parameter_proposal(current_params);
    double alpha_prop = proposed_params[0];

    // Given alpha, update lambda analytically using MLE
    double sum_x_alpha = 0.0;
    arma::vec x = cluster_data.col(0);

    for (int i = 0; i < n; ++i) {
      sum_x_alpha += std::pow(x[i], alpha_prop);
    }

    double lambda_prop = std::pow(sum_x_alpha / n, 1.0 / alpha_prop);
    proposed_params[1] = lambda_prop;

    // Calculate proposed log likelihood
    double proposed_log_lik = 0.0;
    for (int i = 0; i < n; ++i) {
      proposed_log_lik += log_likelihood(cluster_data.row(i).t(), proposed_params);
    }

    double proposed_log_prior = log_prior_density(proposed_params);

    // Calculate acceptance ratio
    double log_ratio = (proposed_log_lik + proposed_log_prior) -
      (current_log_lik + current_log_prior);

    double accept_prob = std::min(1.0, std::exp(log_ratio));

    // Accept or reject
    if (R::runif(0, 1) < accept_prob) {
      current_params = proposed_params;
      current_log_lik = proposed_log_lik;
      current_log_prior = proposed_log_prior;
      accept_count++;
    }
  }

  return current_params;
}

// Update hyperparameters (for hierarchical models)
void WeibullMixing::update_hyperparameters(const std::vector<arma::vec>& all_params) {
  int K = all_params.size();
  if (K == 0) return;

  // Find maximum alpha across all clusters
  double max_alpha = 0.0;
  double sum_inv_lambda = 0.0;

  for (const auto& params : all_params) {
    max_alpha = std::max(max_alpha, params[0]);
    if (params[1] > 1e-10) {
      sum_inv_lambda += 1.0 / params[1];
    }
  }

  // Update phi using Pareto posterior
  double U = R::runif(0, 1);
  phi = qpareto(U, max_alpha, K * hyper_a1 + hyper_a2);

  // Update beta0 using Gamma posterior
  double post_shape = hyper_b1 + K * alpha0;
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
