// src/ExponentialExports.cpp
#include "../inst/include/ExponentialDistribution.h"
#include "../inst/include/RcppConversions.h"

//' @title Draw from an Exponential distribution prior (C++)
//' @description C++ implementation for drawing from the prior distribution of an
//'   Exponential/Gamma model.
//' @param priorParams A numeric vector of prior parameters (alpha0, beta0).
//' @param n The number of samples to draw.
//' @return A list containing the sampled rate parameters (lambda).
//' @export
 // [[Rcpp::export]]
 Rcpp::List exponential_prior_draw_cpp(Rcpp::NumericVector priorParams, int n = 1) {
   dp::ExponentialMixingDistribution md(priorParams);
   return md.priorDraw(n);
 }

//' @title Calculate Exponential log-likelihood (C++)
//' @description C++ implementation for calculating exponential log-likelihood.
//' @param x A numeric vector of data points.
//' @param lambda The rate parameter.
//' @return A numeric vector of log-likelihood values.
//' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector exponential_log_likelihood_cpp(Rcpp::NumericVector x, double lambda) {
   if (lambda <= 0) {
     return Rcpp::NumericVector(x.size(), -std::numeric_limits<double>::infinity());
   }

   const int n = x.size();
   Rcpp::NumericVector result(n);

   for (int i = 0; i < n; ++i) {
     if (x[i] < 0) {
       result[i] = -std::numeric_limits<double>::infinity();
     } else {
       result[i] = std::log(lambda) - lambda * x[i];
     }
   }

   return result;
 }

//' @title Draw from an Exponential distribution posterior (C++)
//' @description C++ implementation for drawing from the posterior distribution of an
//'   Exponential/Gamma model.
//' @param priorParams A numeric vector of prior parameters.
//' @param x A numeric matrix of data points.
//' @param n The number of samples to draw.
//' @return A list containing the sampled rate parameters (lambda).
//' @export
 // [[Rcpp::export]]
 Rcpp::List exponential_posterior_draw_cpp(Rcpp::NumericVector priorParams,
                                           Rcpp::NumericMatrix x,
                                           int n = 1) {
   dp::ExponentialMixingDistribution md(priorParams);
   arma::mat x_arma = Rcpp::as<arma::mat>(x);
   return md.posteriorDraw(x_arma, n);
 }

//' @title Calculate Exponential posterior parameters (C++)
//' @description C++ implementation for calculating posterior parameters for an
//'   Exponential/Gamma model.
//' @param priorParams A numeric vector of prior parameters.
//' @param x A numeric matrix of data.
//' @return A list with alpha and beta posterior parameters.
//' @export
 // [[Rcpp::export]]
 Rcpp::List exponential_posterior_parameters_cpp(Rcpp::NumericVector priorParams,
                                                 Rcpp::NumericMatrix x) {
   dp::ExponentialMixingDistribution md(priorParams);
   arma::mat x_arma = Rcpp::as<arma::mat>(x);

   Rcpp::NumericMatrix post_params = md.posteriorParameters(x_arma);

   return Rcpp::List::create(
     Rcpp::Named("alpha") = post_params(0, 0),
     Rcpp::Named("beta") = post_params(0, 1)
   );
 }

//' @title Calculate Exponential likelihood (C++)
//' @description C++ implementation for calculating exponential likelihood.
//' @param x A numeric vector of data points.
//' @param lambda The rate parameter.
//' @return A numeric vector of likelihood values.
//' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector exponential_likelihood_cpp(Rcpp::NumericVector x, double lambda) {
   if (lambda <= 0) {
     Rcpp::NumericVector result(x.size(), 1e-300);
     return result;
   }

   const int n = x.size();
   Rcpp::NumericVector result(n);

   // Direct pointer access - KEY OPTIMIZATION
   const double* x_ptr = &x[0];
   double* result_ptr = &result[0];

   // Vectorized computation
   for (int i = 0; i < n; ++i) {
     result_ptr[i] = (x_ptr[i] >= 0) ?
     (lambda * std::exp(-lambda * x_ptr[i])) : 0.0;
   }

   return result;
 }

//' @title Calculate Exponential predictive distribution (C++)
//' @description C++ implementation for calculating the predictive distribution.
//' @param priorParams A numeric vector of prior parameters.
//' @param x A numeric vector of data.
//' @return A numeric vector of predictive probabilities.
//' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector exponential_predictive_cpp(Rcpp::NumericVector priorParams,
                                                Rcpp::NumericVector x) {
   dp::ExponentialMixingDistribution md(priorParams);
   arma::vec x_arma = Rcpp::as<arma::vec>(x);
   return md.predictive(x_arma);
 }

//' @title Update cluster components (C++ conjugate exponential)
//' @description C++ implementation of the cluster component update for conjugate models.
//' @param dpObj A list representing the Dirichlet Process object.
//' @return A list with updated cluster assignments and parameters.
//' @export
 // [[Rcpp::export]]
 Rcpp::List conjugate_exponential_cluster_component_update_cpp(Rcpp::List dpObj) {
   // Create C++ DP object using the new constructor
   dp::ConjugateExponentialDP* dp_cpp = new dp::ConjugateExponentialDP(dpObj);

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

//' @title Update alpha for conjugate exponential DP (C++)
//' @description C++ implementation of the concentration parameter update.
//' @param dpObj A list representing the Dirichlet Process object.
//' @return Updated alpha value.
//' @export
 // [[Rcpp::export]]
 double conjugate_exponential_update_alpha_cpp(Rcpp::List dpObj) {
   // Extract necessary components
   double alpha = dpObj["alpha"];
   int n = dpObj["n"];
   int numberClusters = dpObj["numberClusters"];
   Rcpp::NumericVector alphaPriorParameters = dpObj["alphaPriorParameters"];

   // Perform the update using auxiliary variable method
   double x = R::rbeta(alpha + 1.0, n);

   double pi1 = alphaPriorParameters[0] + numberClusters - 1.0;
   double pi2 = n * (alphaPriorParameters[1] - log(x));
   double pi_ratio = pi1 / (pi1 + pi2);

   double postShape, postRate;
   if (R::runif(0, 1) < pi_ratio) {
     postShape = alphaPriorParameters[0] + numberClusters;
   } else {
     postShape = alphaPriorParameters[0] + numberClusters - 1.0;
   }
   postRate = alphaPriorParameters[1] - log(x);

   double new_alpha = R::rgamma(postShape, 1.0/postRate);

   return new_alpha;
 }

//' @title Update cluster parameters (C++ conjugate exponential)
//' @description C++ implementation of the cluster parameter update for conjugate models.
//' @param dpObj A list representing the Dirichlet Process object.
//' @return A list containing the updated cluster parameters.
//' @export
 // [[Rcpp::export]]
 Rcpp::List conjugate_exponential_cluster_parameter_update_cpp(Rcpp::List dpObj) {
   // Create C++ DP object using the new constructor
   dp::ConjugateExponentialDP* dp_cpp = new dp::ConjugateExponentialDP(dpObj);

   // Perform cluster parameter update
   dp_cpp->clusterParameterUpdate();

   // Extract results
   Rcpp::List result = dp_cpp->clusterParameters;

   // Clean up
   delete dp_cpp;

   return result;
 }
