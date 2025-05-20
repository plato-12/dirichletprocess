#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericVector likelihood_normal_cpp(
    const Rcpp::List& mdObj,
    const Rcpp::NumericVector& x,
    const Rcpp::List& theta) {

  // Check that theta has the expected structure
  if (theta.size() != 2) {
    Rcpp::stop("Expected theta to have 2 elements");
  }

  // Extract parameters from theta with more careful conversion
  double mu = Rcpp::as<double>(VECTOR_ELT(theta, 0));
  double sigma = Rcpp::as<double>(VECTOR_ELT(theta, 1));

  int n = x.size();
  Rcpp::NumericVector result(n);

  // Calculate normal likelihood for each point
  for (int i = 0; i < n; i++) {
    double x_i = x[i];

    // Normal density calculation
    double diff = x_i - mu;
    double standardized = diff / sigma;
    result[i] = (1.0 / (sigma * sqrt(2.0 * M_PI))) *
      exp(-0.5 * standardized * standardized);
  }

  return result;
}
