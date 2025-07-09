// src/MVNExports.cpp
#include "../inst/include/MVNormalDistribution.h"
#include "../inst/include/RcppConversions.h"

// The ensureSymmetric function is already defined as inline in the header, so we don't need to define it here

//' @title Draw from a Multivariate Normal-Wishart prior (C++)
 //' @description C++ implementation for drawing from the prior distribution of a
 //'   Multivariate Normal-Wishart model.
 //' @param priorParams A list containing prior parameters (mu0, kappa0, Lambda, nu).
 //' @param n The number of samples to draw.
 //' @return A list containing the sampled parameters (mu and sig).
 //' @export
 // [[Rcpp::export]]
 Rcpp::List mvnormal_prior_draw_cpp(Rcpp::List priorParams, int n = 1) {
   return dp::MVNormalMixingDistribution::priorDrawStatic(priorParams, n);
 }

//' @title Draw from a Multivariate Normal-Wishart posterior (C++)
 //' @description C++ implementation for drawing from the posterior distribution of a
 //'   Multivariate Normal-Wishart model.
 //' @param priorParams A list containing prior parameters.
 //' @param x A numeric matrix of data points.
 //' @param n The number of samples to draw.
 //' @return A list containing the sampled parameters (mu and sig).
 //' @export
 // [[Rcpp::export]]
 Rcpp::List mvnormal_posterior_draw_cpp(Rcpp::List priorParams,
                                        Rcpp::NumericMatrix x,
                                        int n = 1) {
   arma::mat x_arma = Rcpp::as<arma::mat>(x);
   return dp::MVNormalMixingDistribution::posteriorDrawStatic(priorParams, x_arma, n);
 }

//' @title Calculate MVNormal posterior parameters (C++)
 //' @description C++ implementation for calculating posterior parameters for a
 //'   Multivariate Normal-Wishart model.
 //' @param priorParams A list containing prior parameters.
 //' @param x A numeric matrix of data.
 //' @return A list of posterior parameters.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List mvnormal_posterior_parameters_cpp(Rcpp::List priorParams,
                                              Rcpp::NumericMatrix x) {
   dp::MVNormalMixingDistribution md(priorParams);
   arma::mat x_arma = Rcpp::as<arma::mat>(x);
   return md.posteriorParameters(x_arma);
 }

//' @title Calculate MVNormal predictive distribution (C++)
 //' @description C++ implementation for calculating the predictive distribution.
 //' @param priorParams A list containing prior parameters.
 //' @param x A numeric matrix of data.
 //' @return A numeric vector of predictive probabilities.
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector mvnormal_predictive_cpp(Rcpp::List priorParams,
                                             Rcpp::NumericMatrix x) {
   dp::MVNormalMixingDistribution md(priorParams);
   arma::mat x_arma = Rcpp::as<arma::mat>(x);
   return md.predictive(x_arma);
 }

//' @title Calculate MVNormal likelihood (C++)
 //' @description C++ implementation for calculating multivariate normal likelihood.
 //' @param x A numeric matrix of data points.
 //' @param mu Mean vector.
 //' @param sigma Covariance matrix.
 //' @return A numeric vector of likelihood values.
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector mvnormal_likelihood_cpp(Rcpp::NumericMatrix x,
                                             Rcpp::NumericVector mu,
                                             Rcpp::NumericMatrix sigma) {
   arma::mat x_arma = Rcpp::as<arma::mat>(x);
   arma::vec mu_arma = Rcpp::as<arma::vec>(mu);
   arma::mat sigma_arma = Rcpp::as<arma::mat>(sigma);

   // For this export function, sigma is expected to be a covariance matrix
   // (to match mvtnorm::dmvnorm behavior)
   int n = x_arma.n_rows;
   int d = x_arma.n_cols;
   Rcpp::NumericVector result(n);

   // Ensure sigma is symmetric
   sigma_arma = dp::ensureSymmetric(sigma_arma);

   double log_det_val;
   double sign;
   arma::log_det(log_det_val, sign, sigma_arma);

   if (sign <= 0) {
     result.fill(1e-300);
     return result;
   }

   arma::mat sigma_inv;
   try {
     sigma_inv = arma::inv_sympd(sigma_arma);
   } catch(...) {
     result.fill(1e-300);
     return result;
   }

   double log_const = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det_val;

   for (int i = 0; i < n; i++) {
     arma::vec x_centered = x_arma.row(i).t() - mu_arma;
     double quad_form = arma::as_scalar(x_centered.t() * sigma_inv * x_centered);
     result[i] = std::exp(log_const - 0.5 * quad_form);
   }

   // Ensure result is a plain numeric vector without extra attributes
   result.attr("dim") = R_NilValue;
   return result;
 }

//' @title Update cluster components for MVNormal (C++ conjugate)
 //' @description C++ implementation of the cluster component update for MVNormal conjugate models.
 //' @param dpObj A list representing the Dirichlet Process object.
 //' @return A list with updated cluster assignments and parameters.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List conjugate_mvnormal_cluster_component_update_cpp(Rcpp::List dpObj) {
   // Extract necessary components
   arma::mat data = Rcpp::as<arma::mat>(dpObj["data"]);
   arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dpObj["clusterLabels"]);
   arma::uvec pointsPerCluster = Rcpp::as<arma::uvec>(dpObj["pointsPerCluster"]);
   int numberClusters = dpObj["numberClusters"];
   double alpha = dpObj["alpha"];
   Rcpp::List mixingDistribution = dpObj["mixingDistribution"];
   Rcpp::List priorParams = mixingDistribution["priorParameters"];
   Rcpp::List clusterParameters = dpObj["clusterParameters"];
   Rcpp::NumericVector predictiveArray = dpObj["predictiveArray"];

   // Create C++ DP object
   dp::ConjugateMVNormalDP* dp_cpp = new dp::ConjugateMVNormalDP();
   dp_cpp->data = data;
   dp_cpp->n = data.n_rows;
   dp_cpp->alpha = alpha;
   dp_cpp->clusterLabels = clusterLabels;
   dp_cpp->pointsPerCluster = pointsPerCluster;
   dp_cpp->numberClusters = numberClusters;
   dp_cpp->clusterParameters = clusterParameters;
   dp_cpp->predictiveArray = Rcpp::as<arma::vec>(predictiveArray);
   dp_cpp->mixingDistribution = new dp::MVNormalMixingDistribution(priorParams);

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

//' @title Update cluster parameters for MVNormal (C++ conjugate)
 //' @description C++ implementation of the cluster parameter update for MVNormal conjugate models.
 //' @param dpObj A list representing the Dirichlet Process object.
 //' @return A list containing the updated cluster parameters.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List conjugate_mvnormal_cluster_parameter_update_cpp(Rcpp::List dpObj) {
   // Extract necessary components
   arma::mat data = Rcpp::as<arma::mat>(dpObj["data"]);
   arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dpObj["clusterLabels"]);
   int numberClusters = dpObj["numberClusters"];
   Rcpp::List mixingDistribution = dpObj["mixingDistribution"];
   Rcpp::List priorParams = mixingDistribution["priorParameters"];
   Rcpp::List clusterParameters = dpObj["clusterParameters"];

   // Create C++ DP object
   dp::ConjugateMVNormalDP* dp_cpp = new dp::ConjugateMVNormalDP();
   dp_cpp->data = data;
   dp_cpp->n = data.n_rows;
   dp_cpp->clusterLabels = clusterLabels;
   dp_cpp->numberClusters = numberClusters;
   dp_cpp->clusterParameters = clusterParameters;
   dp_cpp->mixingDistribution = new dp::MVNormalMixingDistribution(priorParams);

   // Perform cluster parameter update
   dp_cpp->clusterParameterUpdate();

   // Extract results
   Rcpp::List result = dp_cpp->clusterParameters;

   // Clean up
   delete dp_cpp;

   return result;
 }
