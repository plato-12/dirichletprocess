#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Normal distribution likelihood calculation
// [[Rcpp::export]]
Rcpp::NumericVector normal_likelihood_cpp(const Rcpp::NumericVector& x,
                                          double mu,
                                          double sigma) {
  int n = x.size();
  Rcpp::NumericVector result(n);

  double log_const = -0.5 * std::log(2.0 * M_PI) - std::log(sigma);

  for (int i = 0; i < n; i++) {
    double z = (x[i] - mu) / sigma;
    result[i] = std::exp(log_const - 0.5 * z * z);
  }

  return result;
}

// Generic likelihood dispatcher
// [[Rcpp::export]]
Rcpp::NumericVector likelihood_cpp(const Rcpp::List& mdObj,
                                   const Rcpp::NumericVector& x,
                                   const Rcpp::List& theta) {
  // Extract distribution type from mdObj
  std::string dist_type = Rcpp::as<std::string>(mdObj["distribution"]);

  if (dist_type == "normal") {
    // Extract parameters - handle 3D arrays properly
    Rcpp::NumericVector mu_array = theta[0];
    Rcpp::NumericVector sigma_array = theta[1];

    // Extract the first element from each array
    double mu = mu_array[0];
    double sigma = sigma_array[0];

    return normal_likelihood_cpp(x, mu, sigma);
  } else {
    Rcpp::stop("Distribution type not implemented in C++: " + dist_type);
  }
}

// Multivariate normal distribution likelihood calculation
// [[Rcpp::export]]
arma::vec mvnormal_likelihood_cpp(const arma::mat& x,
                                  const arma::rowvec& mu,
                                  const arma::mat& sigma) {
  int n = x.n_rows;
  arma::vec result(n);

  // Calculate determinant and inverse once
  double log_det;
  double sign;
  arma::log_det(log_det, sign, sigma);
  arma::mat sigma_inv = arma::inv_sympd(sigma);

  double d = x.n_cols;
  double log_const = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det;

  for (int i = 0; i < n; i++) {
    arma::rowvec x_centered = x.row(i) - mu;
    double quad_form = arma::as_scalar(x_centered * sigma_inv * x_centered.t());
    result(i) = std::exp(log_const - 0.5 * quad_form);
  }

  return result;
}
