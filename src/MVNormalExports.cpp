#include <RcppArmadillo.h>
#include "../inst/include/mcmc_runner.h"

// Additional export for testing MVNormal likelihood
// [[Rcpp::export]]
arma::vec mvnormal_log_likelihood_cpp(arma::mat x, arma::vec mu, arma::mat Sigma) {
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec log_lik(n);

  // Ensure Sigma is symmetric
  Sigma = 0.5 * (Sigma + Sigma.t());

  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, Sigma);

  if (sign <= 0) {
    log_lik.fill(-std::numeric_limits<double>::infinity());
    return log_lik;
  }

  arma::mat Sigma_inv;
  try {
    Sigma_inv = arma::inv_sympd(Sigma);
  } catch (...) {
    log_lik.fill(-std::numeric_limits<double>::infinity());
    return log_lik;
  }

  double log_const = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det_val;

  for (int i = 0; i < n; ++i) {
    arma::vec x_centered = x.row(i).t() - mu;
    double quad_form = arma::as_scalar(x_centered.t() * Sigma_inv * x_centered);
    log_lik(i) = log_const - 0.5 * quad_form;
  }

  return log_lik;
}
