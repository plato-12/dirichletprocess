#include <RcppArmadillo.h>
#include "../inst/include/normal_fixed_variance_mixing.h"

using namespace dirichletprocess;

// [[Rcpp::export]]
Rcpp::NumericVector cpp_normal_fixed_variance_prior_draw(double mu0, double sigma0,
                                                         double sigma, int n) {
  NormalFixedVarianceMixing nfv(mu0, sigma0, sigma);
  Rcpp::NumericVector result(n);

  for (int i = 0; i < n; ++i) {
    arma::vec params = nfv.prior_draw();
    result[i] = params[0];
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_normal_fixed_variance_posterior_draw(arma::mat data, double mu0,
                                                             double sigma0, double sigma, int n) {
  NormalFixedVarianceMixing nfv(mu0, sigma0, sigma);
  Rcpp::NumericVector result(n);

  for (int i = 0; i < n; ++i) {
    arma::vec params = nfv.posterior_draw(data, arma::vec());
    result[i] = params[0];
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_normal_fixed_variance_likelihood(arma::vec x, double mu, double sigma) {
  NormalFixedVarianceMixing nfv(0.0, 1.0, sigma);
  arma::vec params(1);
  params[0] = mu;

  Rcpp::NumericVector result(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    result[i] = std::exp(nfv.log_likelihood(x.row(i).t(), params));
  }

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_normal_fixed_variance_posterior_parameters(arma::mat data,
                                                                   double mu0, double sigma0,
                                                                   double sigma) {
  NormalFixedVarianceMixing nfv(mu0, sigma0, sigma);
  arma::vec params = nfv.posterior_parameters(data);

  return Rcpp::NumericVector::create(params[0], params[1]);
}
