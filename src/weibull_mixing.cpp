// src/weibull_mixing.cpp
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
  if (x < 0 || alpha <= 0 || lambda <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Special case for x = 0
  if (x == 0) {
    // For Weibull, when x=0, the density is 0 unless alpha=1 (exponential case)
    return (alpha == 1.0) ? std::log(1.0 / lambda) : -std::numeric_limits<double>::infinity();
  }

  // Weibull PDF: f(x|α,λ) = (α/λ) * (x/λ)^(α-1) * exp(-(x/λ)^α)
  // Log PDF: log(α) - log(λ) + (α-1)*log(x) - (α-1)*log(λ) - (x/λ)^α
  //        = log(α) - α*log(λ) + (α-1)*log(x) - (x/λ)^α
  double log_lik = std::log(alpha) - std::log(lambda) +
    (alpha - 1.0) * std::log(x) -
    std::pow(x, alpha) / lambda;

  return log_lik;
}

// Prior draw
arma::vec WeibullMixing::prior_draw() const {
  arma::vec params(2);

  // alpha ~ Uniform(0, phi)
  params[0] = R::runif(0, phi);

  // lambda = 1/Gamma(alpha0, beta0) where beta0 is the rate parameter
  // R::rgamma uses shape and scale, so we need scale = 1/rate = 1/beta0
  double gamma_draw = R::rgamma(alpha0, 1.0 / beta0);
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

  // Propose new alpha using absolute value to ensure positivity
  proposed[0] = std::abs(current_params[0] + mh_step_alpha * R::rnorm(0.0, 1.7));

  // Bound check
  if (proposed[0] > phi) {
    proposed[0] = phi * R::runif(0.5, 1.0);
  }

  // Lambda remains the same for this proposal (will be updated in posterior_draw)
  return proposed;
}

// Posterior draw using Metropolis-Hastings
arma::vec WeibullMixing::posterior_draw(const arma::mat& cluster_data,
                                        const arma::vec& prior_params) const {
  int n_data = cluster_data.n_rows;

  // Handle empty cluster
  if (n_data == 0) {
    return prior_draw();
  }

  arma::vec x_vec = cluster_data.col(0);

  // Initialize with prior draw
  arma::vec current_params = prior_draw();
  double alpha_current = current_params[0];
  double lambda_current = current_params[1];

  // Pre-compute sum(x^alpha) for current alpha
  double sum_x_alpha_current = 0.0;
  for (int i = 0; i < n_data; i++) {
    if (x_vec[i] > 0) {
      sum_x_alpha_current += std::pow(x_vec[i], alpha_current);
    }
  }

  // Update lambda given alpha (Gibbs step)
  double shape_post = alpha0 + n_data;
  double rate_post_current = sum_x_alpha_current + beta0;
  lambda_current = 1.0 / R::rgamma(shape_post, 1.0 / rate_post_current);

  // Compute initial log likelihood and prior
  double current_log_lik = 0.0;
  for (int i = 0; i < n_data; i++) {
    current_log_lik += log_likelihood(x_vec.row(i).t(),
                                      arma::vec({alpha_current, lambda_current}));
  }
  double current_log_prior = log_prior_density(arma::vec({alpha_current, lambda_current}));

  // MH sampling for alpha
  int accept_count = 0;
  double adaptive_mh_step = mh_step_alpha;

  for (int iter = 0; iter < mh_draws; iter++) {
    // Propose new alpha
    double alpha_prop = std::abs(alpha_current + adaptive_mh_step * R::rnorm(0.0, 1.7));

    // Bound check
    if (alpha_prop > phi) {
      alpha_prop = phi * R::runif(0.5, 1.0);
    }

    // Compute sum(x^alpha) for proposed alpha
    double sum_x_alpha_prop = 0.0;
    bool valid_sum = true;
    for (int i = 0; i < n_data; i++) {
      if (x_vec[i] > 0) {
        double x_alpha = std::pow(x_vec[i], alpha_prop);
        if (std::isfinite(x_alpha)) {
          sum_x_alpha_prop += x_alpha;
        } else {
          valid_sum = false;
          break;
        }
      }
    }

    if (!valid_sum || sum_x_alpha_prop <= 0) {
      continue; // Reject this proposal
    }

    // Sample lambda given proposed alpha
    double rate_post_prop = sum_x_alpha_prop + beta0;
    double lambda_prop = 1.0 / R::rgamma(shape_post, 1.0 / rate_post_prop);

    // Compute proposed log likelihood
    double proposed_log_lik = 0.0;
    for (int i = 0; i < n_data; i++) {
      proposed_log_lik += log_likelihood(x_vec.row(i).t(),
                                         arma::vec({alpha_prop, lambda_prop}));
    }

    // Compute log priors
    double proposed_log_prior = log_prior_density(arma::vec({alpha_prop, lambda_prop}));

    // MH acceptance ratio
    double log_ratio = (proposed_log_lik + proposed_log_prior) -
      (current_log_lik + current_log_prior);

    double accept_prob = std::min(1.0, std::exp(log_ratio));
    if (!std::isfinite(accept_prob)) {
      accept_prob = 0.0;
    }

    // Accept/reject
    if (R::runif(0, 1) < accept_prob) {
      alpha_current = alpha_prop;
      lambda_current = lambda_prop;
      current_log_lik = proposed_log_lik;
      current_log_prior = proposed_log_prior;
      accept_count++;
    }

    // Adaptive step sizing every 50 iterations
    if (iter > 0 && iter % 50 == 0) {
      double recent_accept_rate = (double)accept_count / 50.0;

      if (recent_accept_rate < 0.15) {
        adaptive_mh_step *= 0.7;  // Decrease step size
      } else if (recent_accept_rate > 0.5) {
        adaptive_mh_step *= 1.3;  // Increase step size
      }

      // Reset counter
      accept_count = 0;

      // Keep step size in reasonable bounds
      adaptive_mh_step = std::max(0.01, std::min(5.0, adaptive_mh_step));
    }
  }

  return arma::vec({alpha_current, lambda_current});
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
