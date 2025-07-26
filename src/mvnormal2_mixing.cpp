#include "mvnormal2_mixing.h"
#include <RcppArmadillo.h>
#include <cmath>

namespace dirichletprocess {

MVNormal2Mixing::MVNormal2Mixing(const arma::mat& mu0, const arma::mat& sigma0,
                                 const arma::mat& phi0, double nu0)
  : mu0(mu0), sigma0(sigma0), phi0(phi0), nu0(nu0) {
  
  // Extract dimension from mu0
  if (mu0.n_rows == 1) {
    d = mu0.n_cols;
  } else if (mu0.n_cols == 1) {
    d = mu0.n_rows;
  } else {
    Rcpp::stop("mu0 must be a row or column vector");
  }

  // Validate inputs
  if (nu0 <= d - 1) {
    Rcpp::stop("nu0 must be greater than dimension - 1");
  }
  if (phi0.n_rows != static_cast<arma::uword>(d) || phi0.n_cols != static_cast<arma::uword>(d)) {
    Rcpp::stop("phi0 must be a d x d matrix");
  }
  if (sigma0.n_rows != static_cast<arma::uword>(d) || sigma0.n_cols != static_cast<arma::uword>(d)) {
    Rcpp::stop("sigma0 must be a d x d matrix");
  }

  // Ensure phi0 and sigma0 are symmetric
  this->phi0 = ensureSymmetric(phi0);
  this->sigma0 = ensureSymmetric(sigma0);
}

double MVNormal2Mixing::log_likelihood(const arma::vec& data_point,
                                       const arma::vec& params) const {
  // Extract mean and covariance matrix from flattened params
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

  // Compute quadratic form
  arma::mat Sigma_inv;
  try {
    Sigma_inv = arma::inv_sympd(Sigma);
  } catch (...) {
    return -std::numeric_limits<double>::infinity();
  }

  double quad_form = arma::as_scalar(x_centered.t() * Sigma_inv * x_centered);
  double log_lik = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det_val - 0.5 * quad_form;

  return log_lik;
}

arma::vec MVNormal2Mixing::posterior_draw(const arma::mat& cluster_data,
                                          const arma::vec& prior_params) const {
  int n = cluster_data.n_rows;

  // Handle empty cluster case
  if (n == 0) {
    return prior_draw();
  }

  // MVNormal2 semi-conjugate posterior sampling
  // This is more complex than the fully conjugate case
  
  // For semi-conjugate case, we need to use Gibbs sampling
  // 1. Sample Sigma from Inverse-Wishart given mu and data
  // 2. Sample mu from multivariate normal given Sigma and data
  
  arma::vec x_bar = arma::mean(cluster_data, 0).t();
  
  // Update degrees of freedom
  double nu_n = nu0 + n;
  
  // Compute scatter matrix
  arma::mat S = arma::zeros(d, d);
  for (int i = 0; i < n; ++i) {
    arma::vec xi = cluster_data.row(i).t();
    S += (xi - x_bar) * (xi - x_bar).t();
  }
  
  // Update scale matrix for Inverse-Wishart
  arma::vec mu0_vec;
  if (mu0.n_rows == 1) {
    mu0_vec = mu0.row(0).t();
  } else {
    mu0_vec = mu0.col(0);
  }
  arma::mat phi_n = phi0 + S + (n * sigma0 * arma::inv(sigma0 + n * arma::eye(d, d))) * 
                    (x_bar - mu0_vec) * (x_bar - mu0_vec).t();
  phi_n = ensureSymmetric(phi_n);

  // Sample covariance matrix from Inverse-Wishart
  arma::mat phi_n_inv;
  try {
    phi_n_inv = arma::inv_sympd(phi_n);
  } catch (...) {
    phi_n_inv = arma::pinv(phi_n);
  }
  
  arma::mat Sigma_draw = arma::iwishrnd(phi_n_inv, nu_n);
  Sigma_draw = ensureSymmetric(Sigma_draw);

  // Sample mean from multivariate normal given Sigma
  arma::mat sigma_n_inv = arma::inv_sympd(sigma0) + n * arma::inv_sympd(Sigma_draw);
  arma::mat sigma_n = arma::inv_sympd(sigma_n_inv);
  arma::vec mu_n = sigma_n * (arma::inv_sympd(sigma0) * mu0_vec + 
                              n * arma::inv_sympd(Sigma_draw) * x_bar);

  arma::vec mu_draw = arma::mvnrnd(mu_n, sigma_n);

  return flatten_params(mu_draw, Sigma_draw);
}

arma::vec MVNormal2Mixing::prior_draw() const {
  // Draw from prior: mu ~ N(mu0, sigma0), Sigma ~ IW(phi0, nu0)
  
  arma::vec mu0_vec;
  if (mu0.n_rows == 1) {
    mu0_vec = mu0.row(0).t();
  } else {
    mu0_vec = mu0.col(0);
  }
  
  // Sample covariance matrix from Inverse-Wishart
  arma::mat phi0_inv;
  try {
    phi0_inv = arma::inv_sympd(phi0);
  } catch (...) {
    phi0_inv = arma::pinv(phi0);
  }
  
  arma::mat Sigma_draw = arma::iwishrnd(phi0_inv, nu0);
  Sigma_draw = ensureSymmetric(Sigma_draw);

  // Sample mean from multivariate normal
  arma::vec mu_draw = arma::mvnrnd(mu0_vec, sigma0);

  return flatten_params(mu_draw, Sigma_draw);
}

int MVNormal2Mixing::param_dim() const {
  // d parameters for mu + d*(d+1)/2 parameters for symmetric Sigma
  return d + d * (d + 1) / 2;
}

arma::vec MVNormal2Mixing::flatten_params(const arma::vec& mu, const arma::mat& Sigma) const {
  arma::vec params(param_dim());
  
  // First d elements are mu
  params.subvec(0, d - 1) = mu;
  
  // Remaining elements are upper triangle of Sigma (including diagonal)
  int idx = d;
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      params(idx++) = Sigma(i, j);
    }
  }
  
  return params;
}

void MVNormal2Mixing::unflatten_params(const arma::vec& params, arma::vec& mu, arma::mat& Sigma) const {
  // Extract mu
  mu = params.subvec(0, d - 1);
  
  // Extract Sigma from upper triangle
  Sigma = arma::zeros(d, d);
  int idx = d;
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      Sigma(i, j) = params(idx);
      if (i != j) {
        Sigma(j, i) = params(idx);  // Make symmetric
      }
      idx++;
    }
  }
}

} // namespace dirichletprocess