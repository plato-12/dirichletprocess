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

  // Weibull PDF: f(x|α,λ) = (1/λ) * α * (x)^(α-1) * exp(-(1/λ) * x^α)
  // Note: This matches the R parameterization where lambda appears as 1/lambda
  double log_lik = -std::log(lambda) + std::log(alpha) +
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

  // Propose new alpha (matching R's abs() approach)
  double alpha_current = current_params[0];
  double alpha_proposal = std::abs(alpha_current + mh_step_alpha * R::rnorm(0, 1.7));

  // Ensure alpha stays within bounds [0, phi]
  if (alpha_proposal > phi) {
    alpha_proposal = phi;
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

  // Run Metropolis-Hastings with Gibbs update for lambda
  for (int iter = 0; iter < mh_draws; ++iter) {
    double alpha_current = current_params[0];
    double lambda_current = current_params[1];

    // Propose new alpha
    double alpha_prop = std::abs(alpha_current + mh_step_alpha * R::rnorm(0, 1.7));
    if (alpha_prop > phi) {
      alpha_prop = phi;
    }

    // Given current alpha, update lambda using conjugate posterior
    // This matches the R implementation in MetropolisHastings.weibull
    double sum_x_alpha_current = 0.0;
    for (int i = 0; i < n; ++i) {
      double xi = cluster_data(i, 0);
      if (xi > 0) {
        sum_x_alpha_current += std::pow(xi, alpha_current);
      }
    }

    // Sample new lambda from inverse gamma (conjugate posterior)
    // lambda | data, alpha ~ InvGamma(n + alpha0, sum(x^alpha) + beta0)
    double shape_post = n + alpha0;
    double rate_post = sum_x_alpha_current + beta0;
    double gamma_sample = R::rgamma(shape_post, 1.0 / rate_post);
    lambda_current = 1.0 / std::max(1e-10, gamma_sample);

    // Now compute sum for proposed alpha
    double sum_x_alpha_prop = 0.0;
    for (int i = 0; i < n; ++i) {
      double xi = cluster_data(i, 0);
      if (xi > 0) {
        sum_x_alpha_prop += std::pow(xi, alpha_prop);
      }
    }

    // Sample lambda for proposed alpha
    double gamma_sample_prop = R::rgamma(shape_post, 1.0 / (sum_x_alpha_prop + beta0));
    double lambda_prop = 1.0 / std::max(1e-10, gamma_sample_prop);

    // Compute log likelihoods
    double log_lik_current = 0.0;
    double log_lik_prop = 0.0;

    arma::vec params_current = {alpha_current, lambda_current};
    arma::vec params_prop = {alpha_prop, lambda_prop};

    for (int i = 0; i < n; ++i) {
      arma::vec xi = cluster_data.row(i).t();
      log_lik_current += log_likelihood(xi, params_current);
      log_lik_prop += log_likelihood(xi, params_prop);
    }

    // Compute log priors (only for alpha since lambda is integrated out)
    double log_prior_alpha_current = (alpha_current > 0 && alpha_current <= phi) ?
    -std::log(phi) : -std::numeric_limits<double>::infinity();
    double log_prior_alpha_prop = (alpha_prop > 0 && alpha_prop <= phi) ?
    -std::log(phi) : -std::numeric_limits<double>::infinity();

    // Accept/reject for alpha
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
  // The R code uses: rpareto(n, max(clusterParameters[[1]], hyperPriorParameters[1]),
  //                          hyperPriorParameters[2] + numClusters)
  double xm = std::max(max_alpha, hyper_a1);
  double shape = hyper_a2 + K;
  double U = R::runif(0, 1);
  phi = qpareto(U, xm, shape);

  // Update beta0 using Gamma posterior
  // The R code uses: rgamma(n, hyperPriorParameters[3] + 2 * numClusters,
  //                        hyperPriorParameters[4] + sum(1/clusterParameters[[2]]))
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
