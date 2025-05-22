// src/NormalDistribution.cpp
#include "../inst/include/NormalDistribution.h"
#include "../inst/include/RcppConversions.h"

namespace dp {

// NormalMixingDistribution implementation
NormalMixingDistribution::NormalMixingDistribution(const Rcpp::NumericVector& priorParams) {
  distribution = "normal";
  conjugate = true;
  priorParameters = priorParams;
}

NormalMixingDistribution::~NormalMixingDistribution() {
  // Destructor
}

Rcpp::NumericVector NormalMixingDistribution::likelihood(const arma::vec& x, const Rcpp::List& theta) const {
  // Extract parameters from theta
  Rcpp::NumericVector mu_array = theta[0];
  Rcpp::NumericVector sigma_array = theta[1];

  int n_clusters = mu_array.size();
  int n_data = x.n_elem;

  Rcpp::NumericVector result(n_data * n_clusters);

  for (int k = 0; k < n_clusters; k++) {
    double mu = mu_array[k];
    double sigma = sigma_array[k];

    for (int i = 0; i < n_data; i++) {
      result[i * n_clusters + k] = R::dnorm(x[i], mu, sigma, false);
    }
  }

  return result;
}

Rcpp::List NormalMixingDistribution::priorDraw(int n) const {
  Rcpp::NumericVector mu(n);
  Rcpp::NumericVector sigma(n);

  for (int i = 0; i < n; i++) {
    double lambda = R::rgamma(priorParameters[2], 1.0/priorParameters[3]);
    mu[i] = R::rnorm(priorParameters[0], 1.0/sqrt(priorParameters[1] * lambda));
    sigma[i] = sqrt(1.0/lambda);
  }

  // Convert to 3D arrays
  Rcpp::NumericVector mu_arr(n);
  Rcpp::NumericVector sigma_arr(n);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);
  sigma_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    mu_arr[i] = mu[i];
    sigma_arr[i] = sigma[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sigma") = sigma_arr
  );
}

Rcpp::List NormalMixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  Rcpp::NumericMatrix postParams = posteriorParameters(x);

  Rcpp::NumericVector mu(n);
  Rcpp::NumericVector sigma(n);

  for (int i = 0; i < n; i++) {
    double lambda = R::rgamma(postParams(0, 2), 1.0/postParams(0, 3));
    mu[i] = R::rnorm(postParams(0, 0), 1.0/sqrt(postParams(0, 1) * lambda));
    sigma[i] = sqrt(1.0/lambda);
  }

  // Convert to 3D arrays
  Rcpp::NumericVector mu_arr(n);
  Rcpp::NumericVector sigma_arr(n);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);
  sigma_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    mu_arr[i] = mu[i];
    sigma_arr[i] = sigma[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sigma") = sigma_arr
  );
}

Rcpp::NumericMatrix NormalMixingDistribution::posteriorParameters(const arma::mat& x) const {
  int n_x = x.n_rows;
  double ybar = arma::mean(arma::vectorise(x));

  double mu0 = priorParameters[0];
  double kappa0 = priorParameters[1];
  double alpha0 = priorParameters[2];
  double beta0 = priorParameters[3];

  double mu_n = (kappa0 * mu0 + n_x * ybar) / (kappa0 + n_x);
  double kappa_n = kappa0 + n_x;
  double alpha_n = alpha0 + n_x / 2.0;
  double beta_n = beta0 + 0.5 * arma::sum(arma::square(arma::vectorise(x) - ybar)) +
    kappa0 * n_x * std::pow(ybar - mu0, 2) / (2.0 * (kappa0 + n_x));

  Rcpp::NumericMatrix result(1, 4);
  result(0, 0) = mu_n;
  result(0, 1) = kappa_n;
  result(0, 2) = alpha_n;
  result(0, 3) = beta_n;

  return result;
}

Rcpp::NumericVector NormalMixingDistribution::predictive(const arma::vec& x) const {
  int n = x.n_elem;
  Rcpp::NumericVector result(n);

  for (int i = 0; i < n; i++) {
    Rcpp::NumericMatrix postParams = posteriorParameters(arma::mat(&x[i], 1, 1));

    double predictive_val = (R::gammafn(postParams(0, 2)) / R::gammafn(priorParameters[2])) *
      (std::pow(priorParameters[3], priorParameters[2]) /
      std::pow(postParams(0, 3), postParams(0, 2))) *
      std::sqrt(priorParameters[1] / postParams(0, 1));

    result[i] = predictive_val;
  }

  return result;
}

// ConjugateNormalDP implementation
ConjugateNormalDP::ConjugateNormalDP() : mixingDistribution(nullptr), numberClusters(0) {
  // Constructor
}

ConjugateNormalDP::~ConjugateNormalDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void ConjugateNormalDP::clusterComponentUpdate() {
  // Implementation of conjugate cluster component update
  // This is a simplified version - full implementation would be more complex
  Rcpp::warning("ConjugateNormalDP::clusterComponentUpdate not fully implemented");
}

void ConjugateNormalDP::clusterParameterUpdate() {
  // Implementation of conjugate cluster parameter update
  Rcpp::warning("ConjugateNormalDP::clusterParameterUpdate not fully implemented");
}

void ConjugateNormalDP::updateAlpha() {
  // Implementation of alpha update
  Rcpp::warning("ConjugateNormalDP::updateAlpha not fully implemented");
}

Rcpp::List ConjugateNormalDP::clusterLabelChange(int i, int newLabel, int currentLabel) {
  // Implementation of cluster label change
  Rcpp::warning("ConjugateNormalDP::clusterLabelChange not fully implemented");
  return Rcpp::List::create();
}

void ConjugateNormalDP::initialisePredictive() {
  // Implementation of predictive initialization
  Rcpp::warning("ConjugateNormalDP::initialisePredictive not fully implemented");
}

} // namespace dp
