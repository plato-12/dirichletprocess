// src/NormalExports.cpp
#include "../inst/include/NormalDistribution.h"
#include "../inst/include/RcppConversions.h"

// [[Rcpp::export]]
Rcpp::List normal_prior_draw_cpp(Rcpp::NumericVector priorParams, int n = 1) {
  return dp::NormalMixingDistribution::priorDrawStatic(priorParams, n);
}

// [[Rcpp::export]]
Rcpp::List normal_posterior_draw_cpp(Rcpp::NumericVector priorParams,
                                     Rcpp::NumericMatrix x,
                                     int n = 1) {
  arma::mat x_arma = Rcpp::as<arma::mat>(x);
  return dp::NormalMixingDistribution::posteriorDrawStatic(priorParams, x_arma, n);
}

// [[Rcpp::export]]
Rcpp::List conjugate_cluster_component_update_cpp(Rcpp::List dpObj) {
  // Extract necessary components from dpObj
  arma::mat data = Rcpp::as<arma::mat>(dpObj["data"]);
  arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dpObj["clusterLabels"]);
  arma::uvec pointsPerCluster = Rcpp::as<arma::uvec>(dpObj["pointsPerCluster"]);
  int numberClusters = dpObj["numberClusters"];
  double alpha = dpObj["alpha"];
  Rcpp::List mixingDistribution = dpObj["mixingDistribution"];
  Rcpp::NumericVector priorParams = mixingDistribution["priorParameters"];
  Rcpp::List clusterParameters = dpObj["clusterParameters"];
  Rcpp::NumericVector predictiveArray = dpObj["predictiveArray"];

  // Create C++ DP object
  dp::ConjugateNormalDP* dp_cpp = new dp::ConjugateNormalDP();
  dp_cpp->data = data;
  dp_cpp->n = data.n_rows;
  dp_cpp->alpha = alpha;
  dp_cpp->clusterLabels = clusterLabels;
  dp_cpp->pointsPerCluster = pointsPerCluster;
  dp_cpp->numberClusters = numberClusters;
  dp_cpp->clusterParameters = clusterParameters;
  dp_cpp->predictiveArray = Rcpp::as<arma::vec>(predictiveArray);
  dp_cpp->mixingDistribution = new dp::NormalMixingDistribution(priorParams);

  // Perform cluster component update
  dp_cpp->clusterComponentUpdate();

  // Extract results
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("clusterLabels") = dp_cpp->clusterLabels,
    Rcpp::Named("pointsPerCluster") = dp_cpp->pointsPerCluster,
    Rcpp::Named("numberClusters") = dp_cpp->numberClusters,
    Rcpp::Named("clusterParameters") = dp_cpp->clusterParameters
  );

  // Clean up
  delete dp_cpp;

  return result;
}

// [[Rcpp::export]]
Rcpp::List conjugate_cluster_parameter_update_cpp(Rcpp::List dpObj) {
  // Extract necessary components from dpObj
  arma::mat data = Rcpp::as<arma::mat>(dpObj["data"]);
  arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dpObj["clusterLabels"]);
  int numberClusters = dpObj["numberClusters"];
  Rcpp::List mixingDistribution = dpObj["mixingDistribution"];
  Rcpp::NumericVector priorParams = mixingDistribution["priorParameters"];
  Rcpp::List clusterParameters = dpObj["clusterParameters"];

  // Create C++ DP object
  dp::ConjugateNormalDP* dp_cpp = new dp::ConjugateNormalDP();
  dp_cpp->data = data;
  dp_cpp->n = data.n_rows;
  dp_cpp->clusterLabels = clusterLabels;
  dp_cpp->numberClusters = numberClusters;
  dp_cpp->clusterParameters = clusterParameters;
  dp_cpp->mixingDistribution = new dp::NormalMixingDistribution(priorParams);

  // Perform cluster parameter update
  dp_cpp->clusterParameterUpdate();

  // Extract results
  Rcpp::List result = dp_cpp->clusterParameters;

  // Clean up
  delete dp_cpp;

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix normal_posterior_parameters_cpp(Rcpp::NumericVector priorParams,
                                                    Rcpp::NumericMatrix x) {
  dp::NormalMixingDistribution md(priorParams);
  arma::mat x_arma = Rcpp::as<arma::mat>(x);
  return md.posteriorParameters(x_arma);
}
