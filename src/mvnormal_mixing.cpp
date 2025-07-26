#include "mvnormal_mixing.h"
#include <RcppArmadillo.h>
#include <cmath>

namespace dirichletprocess {

MVNormalMixing::MVNormalMixing(const arma::vec& mu0, double kappa0,
                               const arma::mat& Lambda, double nu)
  : mu0(mu0), kappa0(kappa0), Lambda(Lambda), nu(nu), d(mu0.n_elem) {

  // Validate inputs
  if (kappa0 <= 0) {
    Rcpp::stop("kappa0 must be positive");
  }
  if (nu <= d - 1) {
    Rcpp::stop("nu must be greater than dimension - 1");
  }
  if (Lambda.n_rows != d || Lambda.n_cols != d) {
    Rcpp::stop("Lambda must be a d x d matrix");
  }

  // Ensure Lambda is symmetric
  this->Lambda = ensureSymmetric(Lambda);
}

double MVNormalMixing::log_likelihood(const arma::vec& data_point,
                                      const arma::vec& params) const {
  // Extract mean and precision matrix from flattened params
  arma::vec mu(d);
  arma::mat Sigma(d, d);
  unflatten_params(params, mu, Sigma);

  // Compute log-likelihood for multivariate normal
  arma::vec x_centered = data_point - mu;

  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, Sigma);

  if (sign <= 0) {
    return -std::numeric_limits<double>::infinity();
  }

  // Using precision parameterization
  double quad_form = arma::as_scalar(x_centered.t() * Sigma * x_centered);
  double log_lik = -0.5 * d * std::log(2.0 * M_PI) + 0.5 * log_det_val - 0.5 * quad_form;

  return log_lik;
}

arma::vec MVNormalMixing::posterior_draw(const arma::mat& cluster_data,
                                         const arma::vec& prior_params) const {
  int n = cluster_data.n_rows;

  // Handle empty cluster case
  if (n == 0) {
    return prior_draw();
  }

  // Compute sample statistics
  arma::vec x_bar = arma::mean(cluster_data, 0).t();

  // Update posterior parameters (Normal-Wishart conjugate update)
  double kappa_n = kappa0 + n;
  arma::vec mu_n = (kappa0 * mu0 + n * x_bar) / kappa_n;
  double nu_n = nu + n;

  // Compute scatter matrix
  arma::mat S = arma::zeros(d, d);
  for (int i = 0; i < n; ++i) {
    arma::vec xi = cluster_data.row(i).t();
    S += (xi - x_bar) * (xi - x_bar).t();
  }

  // Update Lambda_n
  arma::mat Lambda_n = Lambda + S +
    (kappa0 * n / kappa_n) * (x_bar - mu0) * (x_bar - mu0).t();
  Lambda_n = ensureSymmetric(Lambda_n);

  // Draw from posterior Wishart(nu_n, Lambda_n^{-1})
  arma::mat Lambda_n_inv;
  try {
    Lambda_n_inv = arma::inv_sympd(Lambda_n);
  } catch (...) {
    Rcpp::warning("Matrix inversion failed in posterior_draw, using pseudo-inverse");
    Lambda_n_inv = arma::pinv(Lambda_n);
  }

  // Sample precision matrix from Wishart distribution
  arma::mat prec_draw = arma::wishrnd(Lambda_n_inv, nu_n);
  prec_draw = ensureSymmetric(prec_draw);

  // Sample mean from multivariate normal
  arma::mat prec_mu = kappa_n * prec_draw;
  arma::mat cov_mu;
  try {
    cov_mu = arma::inv_sympd(prec_mu);
  } catch (...) {
    cov_mu = arma::pinv(prec_mu);
  }

  arma::vec mu_draw = arma::mvnrnd(mu_n, cov_mu);

  // Return flattened parameters
  return flatten_params(mu_draw, prec_draw);
}

arma::vec MVNormalMixing::prior_draw() const {
  // Draw precision matrix from Wishart(nu, Lambda^{-1})
  arma::mat Lambda_inv;
  try {
    Lambda_inv = arma::inv_sympd(Lambda);
  } catch (...) {
    Rcpp::warning("Matrix inversion failed in prior_draw, using pseudo-inverse");
    Lambda_inv = arma::pinv(Lambda);
  }

  arma::mat prec_draw = arma::wishrnd(Lambda_inv, nu);
  prec_draw = ensureSymmetric(prec_draw);

  // Draw mean from multivariate normal
  arma::mat prec_mu = kappa0 * prec_draw;
  arma::mat cov_mu;
  try {
    cov_mu = arma::inv_sympd(prec_mu);
  } catch (...) {
    cov_mu = arma::pinv(prec_mu);
  }

  arma::vec mu_draw = arma::mvnrnd(mu0, cov_mu);

  return flatten_params(mu_draw, prec_draw);
}

int MVNormalMixing::param_dim() const {
  // Mean vector (d) + precision matrix (d*d)
  return d + d * d;
}

arma::vec MVNormalMixing::flatten_params(const arma::vec& mu,
                                         const arma::mat& Sigma) const {
  arma::vec params(param_dim());

  // First d elements are the mean
  params.subvec(0, d-1) = mu;

  // Remaining elements are the precision matrix (column-major order)
  params.subvec(d, param_dim()-1) = arma::vectorise(Sigma);

  return params;
}

void MVNormalMixing::unflatten_params(const arma::vec& params,
                                      arma::vec& mu, arma::mat& Sigma) const {
  // Extract mean
  mu = params.subvec(0, d-1);

  // Extract precision matrix
  arma::vec sigma_vec = params.subvec(d, param_dim()-1);
  Sigma = arma::reshape(sigma_vec, d, d);
  Sigma = ensureSymmetric(Sigma);
}

double MVNormalMixing::predictive_probability(const arma::vec& data_point) const {
  // Normal-Wishart predictive distribution
  // This is the multivariate Student's t-distribution
  
  arma::vec x = data_point;
  
  // Predictive parameters
  arma::vec mu_pred = mu0;
  double nu_pred = nu - d + 1;
  
  // Scale matrix for predictive distribution
  arma::mat Lambda_inv;
  try {
    Lambda_inv = arma::inv_sympd(Lambda);
  } catch (...) {
    Lambda_inv = arma::pinv(Lambda);
  }
  
  arma::mat Scale = Lambda_inv * (kappa0 + 1) / (kappa0 * nu_pred);
  
  // Compute multivariate t log-density
  arma::vec x_centered = x - mu_pred;
  
  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, Scale);
  
  if (sign <= 0) {
    return 0.0;  // Invalid covariance
  }
  
  arma::mat Scale_inv;
  try {
    Scale_inv = arma::inv_sympd(Scale);
  } catch (...) {
    return 0.0;
  }
  
  double quad_form = arma::as_scalar(x_centered.t() * Scale_inv * x_centered);
  
  // Log probability of multivariate t-distribution
  double log_prob = std::lgamma((nu_pred + d) / 2.0) - std::lgamma(nu_pred / 2.0) - 
                   (d / 2.0) * std::log(nu_pred * M_PI) - 0.5 * log_det_val - 
                   ((nu_pred + d) / 2.0) * std::log(1 + quad_form / nu_pred);
  
  return std::exp(log_prob);
}

} // namespace dirichletprocess
