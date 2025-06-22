// src/exponential_exports.cpp
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "../inst/include/ExponentialDistribution.h"

using namespace dirichletprocess;

// [[Rcpp::export]]
Rcpp::NumericVector exponential_log_likelihood_cpp(
    Rcpp::NumericVector x, double lambda) {

  int n = x.size();
  Rcpp::NumericVector log_lik(n);

  ExponentialMixing exp_dist(1.0, 1.0); // Dummy prior params
  arma::vec params(1);
  params[0] = lambda;

  for (int i = 0; i < n; i++) {
    arma::vec data_point(1);
    data_point[0] = x[i];
    log_lik[i] = exp_dist.log_likelihood(data_point, params);
  }

  return log_lik;
}

// [[Rcpp::export]]
Rcpp::List exponential_posterior_parameters_cpp(
    Rcpp::NumericVector prior_params,
    Rcpp::NumericMatrix data) {

  double alpha0 = prior_params[0];
  double beta0 = prior_params[1];

  ExponentialMixing exp_dist(alpha0, beta0);
  arma::mat data_arma = Rcpp::as<arma::mat>(data);

  arma::vec post_params = exp_dist.posterior_parameters(data_arma);

  return Rcpp::List::create(
    Rcpp::Named("alpha") = post_params[0],
                                      Rcpp::Named("beta") = post_params[1]
  );
}

// [[Rcpp::export]]
double exponential_posterior_draw_cpp(
    Rcpp::NumericVector prior_params,
    Rcpp::NumericMatrix data) {

  double alpha0 = prior_params[0];
  double beta0 = prior_params[1];

  ExponentialMixing exp_dist(alpha0, beta0);
  arma::mat data_arma = Rcpp::as<arma::mat>(data);
  arma::vec dummy_prior(0); // Not used for conjugate case

  arma::vec params = exp_dist.posterior_draw(data_arma, dummy_prior);
  return params[0];
}

// [[Rcpp::export]]
double exponential_prior_draw_cpp(Rcpp::NumericVector prior_params) {
  double alpha0 = prior_params[0];
  double beta0 = prior_params[1];

  ExponentialMixing exp_dist(alpha0, beta0);
  arma::vec params = exp_dist.prior_draw();
  return params[0];
}
