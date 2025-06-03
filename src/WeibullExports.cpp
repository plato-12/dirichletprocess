// src/WeibullExports.cpp
#include <RcppArmadillo.h>
#include "../inst/include/WeibullDistribution.h"
#include "../inst/include/RcppConversions.h"

// [[Rcpp::export]]
Rcpp::List weibull_prior_draw_cpp(const Rcpp::NumericVector& priorParams, int n = 1) {
  dp::WeibullMixingDistribution md(priorParams, Rcpp::NumericVector::create(1.0, 1.0));
  return md.priorDraw(n);
}

// [[Rcpp::export]]
Rcpp::NumericVector weibull_likelihood_cpp(const Rcpp::NumericVector& x, double alpha, double lambda) {
  arma::vec x_arma = Rcpp::as<arma::vec>(x);

  Rcpp::NumericVector alpha_arr(1);
  Rcpp::NumericVector lambda_arr(1);
  alpha_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  lambda_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  alpha_arr[0] = alpha;
  lambda_arr[0] = lambda;

  Rcpp::List theta = Rcpp::List::create(alpha_arr, lambda_arr);

  dp::WeibullMixingDistribution md(Rcpp::NumericVector::create(10.0, 2.0, 4.0),
                                   Rcpp::NumericVector::create(1.0, 1.0));
  return md.likelihood(x_arma, theta);
}

// [[Rcpp::export]]
double weibull_prior_density_cpp(double alpha, const Rcpp::NumericVector& priorParams) {
  dp::WeibullMixingDistribution md(priorParams, Rcpp::NumericVector::create(1.0, 1.0));

  Rcpp::NumericVector alpha_arr(1);
  Rcpp::NumericVector lambda_arr(1);
  alpha_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  lambda_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  alpha_arr[0] = alpha;
  lambda_arr[0] = 1.0; // Lambda doesn't affect uniform prior on alpha

  Rcpp::List theta = Rcpp::List::create(alpha_arr, lambda_arr);
  return md.priorDensity(theta)[0];
}

// [[Rcpp::export]]
Rcpp::List weibull_posterior_draw_cpp(const Rcpp::NumericVector& priorParams,
                                      const Rcpp::NumericVector& mhStepSize,
                                      const Rcpp::NumericMatrix& x,
                                      int n = 1) {
  arma::mat x_arma = Rcpp::as<arma::mat>(x);
  dp::WeibullMixingDistribution md(priorParams, mhStepSize);
  return md.posteriorDraw(x_arma, n);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix weibull_prior_parameters_update_cpp(const Rcpp::NumericVector& priorParams,
                                                        const Rcpp::NumericVector& hyperPriorParams,
                                                        const Rcpp::List& clusterParameters,
                                                        int n = 1) {
  dp::WeibullMixingDistribution md(priorParams, Rcpp::NumericVector::create(1.0, 1.0), hyperPriorParams);
  md.updatePriorParameters(clusterParameters, n);
  return Rcpp::as<Rcpp::NumericMatrix>(md.priorParameters);
}

// [[Rcpp::export]]
Rcpp::List nonconjugate_weibull_cluster_parameter_update_cpp(Rcpp::List dp_list) {
  // Extract necessary components
  arma::mat data = Rcpp::as<arma::mat>(dp_list["data"]);
  arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dp_list["clusterLabels"]);
  int numberClusters = dp_list["numberClusters"];
  Rcpp::List mixingDistribution = dp_list["mixingDistribution"];
  Rcpp::NumericVector priorParams = mixingDistribution["priorParameters"];
  Rcpp::NumericVector mhStepSize = mixingDistribution["mhStepSize"];
  Rcpp::List clusterParameters = dp_list["clusterParameters"];

  // Create C++ DP object
  dp::NonConjugateWeibullDP* dp_cpp = new dp::NonConjugateWeibullDP();
  dp_cpp->data = data;
  dp_cpp->n = data.n_rows;
  dp_cpp->clusterLabels = clusterLabels;
  dp_cpp->numberClusters = numberClusters;
  dp_cpp->clusterParameters = clusterParameters;

  Rcpp::NumericVector hyperPriorParams;
  if (mixingDistribution.containsElementNamed("hyperPriorParameters")) {
    hyperPriorParams = Rcpp::as<Rcpp::NumericVector>(mixingDistribution["hyperPriorParameters"]);
  }

  dp_cpp->mixingDistribution = new dp::WeibullMixingDistribution(priorParams, mhStepSize, hyperPriorParams);

  // Perform cluster parameter update
  dp_cpp->clusterParameterUpdate();

  // Extract results
  Rcpp::List result = dp_cpp->clusterParameters;

  // Clean up
  delete dp_cpp;

  return result;
}

// [[Rcpp::export]]
Rcpp::List nonconjugate_weibull_cluster_component_update_cpp(Rcpp::List dp_list) {
  // Extract necessary components
  arma::mat data = Rcpp::as<arma::mat>(dp_list["data"]);
  arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dp_list["clusterLabels"]);
  arma::uvec pointsPerCluster = Rcpp::as<arma::uvec>(dp_list["pointsPerCluster"]);
  int numberClusters = dp_list["numberClusters"];
  double alpha = dp_list["alpha"];
  Rcpp::List mixingDistribution = dp_list["mixingDistribution"];
  Rcpp::NumericVector priorParams = mixingDistribution["priorParameters"];
  Rcpp::NumericVector mhStepSize = mixingDistribution["mhStepSize"];
  Rcpp::List clusterParameters = dp_list["clusterParameters"];
  Rcpp::NumericVector alphaPriorParameters = dp_list["alphaPriorParameters"];
  int m = dp_list["m"];

  // Create C++ DP object
  dp::NonConjugateWeibullDP* dp_cpp = new dp::NonConjugateWeibullDP();
  dp_cpp->data = data;
  dp_cpp->n = data.n_rows;
  dp_cpp->alpha = alpha;
  dp_cpp->clusterLabels = clusterLabels;
  dp_cpp->pointsPerCluster = pointsPerCluster;
  dp_cpp->numberClusters = numberClusters;
  dp_cpp->clusterParameters = clusterParameters;
  dp_cpp->alphaPriorParameters = alphaPriorParameters;
  dp_cpp->m = m;

  Rcpp::NumericVector hyperPriorParams;
  if (mixingDistribution.containsElementNamed("hyperPriorParameters")) {
    hyperPriorParams = Rcpp::as<Rcpp::NumericVector>(mixingDistribution["hyperPriorParameters"]);
  }

  dp_cpp->mixingDistribution = new dp::WeibullMixingDistribution(priorParams, mhStepSize, hyperPriorParams);

  // Perform cluster component update
  dp_cpp->clusterComponentUpdate();

  // Convert cluster labels back to 1-indexed for R
  arma::uvec clusterLabels_r = dp_cpp->clusterLabels + 1;

  // Extract results
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("clusterLabels") = clusterLabels_r,  // Convert to 1-indexed
    Rcpp::Named("pointsPerCluster") = dp_cpp->pointsPerCluster,
    Rcpp::Named("numberClusters") = dp_cpp->numberClusters,
    Rcpp::Named("clusterParameters") = dp_cpp->clusterParameters
  );

  // Clean up
  delete dp_cpp;

  return result;
}
