// src/BetaExports.cpp
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"

//' @title Beta distribution likelihood (C++)
 //' @description C++ implementation for calculating the likelihood of data under a Beta distribution
 //' @param x A numeric vector of data points
 //' @param mu The mean parameter of the Beta distribution
 //' @param nu The precision parameter of the Beta distribution
 //' @param maxT The upper bound of the Beta distribution
 //' @return A numeric vector of likelihood values
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector beta_likelihood_cpp(const Rcpp::NumericVector& x_data, // Renamed x
                                         double mu_val, double nu_val, double maxT_val) { // Renamed params
   arma::vec x_arma = Rcpp::as<arma::vec>(x_data);
   return dp::BetaMixingDistribution::likelihoodStatic(x_arma, mu_val, nu_val, maxT_val);
 }

//' @title Draw from a Beta distribution prior (C++)
 //' @description C++ implementation for drawing from the prior distribution of a Beta model
 //' @param priorParams A numeric vector of prior parameters (shape and rate for inverse gamma on nu)
 //' @param maxT The upper bound of the Beta distribution
 //' @param n The number of samples to draw
 //' @return A list containing the sampled parameters (mu and nu)
 //' @export
 // [[Rcpp::export]]
 Rcpp::List beta_prior_draw_cpp(const Rcpp::NumericVector& priorParams,
                                double maxT_val, int n_draws = 1) { // Renamed params
   return dp::BetaMixingDistribution::priorDrawStatic(priorParams, maxT_val, n_draws);
 }

//' @title Draw from a Beta distribution posterior (C++)
 //' @description C++ implementation for drawing from the posterior distribution of a Beta model using Metropolis-Hastings
 //' @param priorParams A numeric vector of prior parameters
 //' @param maxT The upper bound of the Beta distribution
 //' @param mhStepSize A numeric vector of step sizes for the MH algorithm (for mu and nu)
 //' @param x A numeric matrix of data points
 //' @param n The number of samples to return (last n samples from MH chain)
 //' @param mhDrawsVal The total number of Metropolis-Hastings iterations
 //' @return A list containing the sampled parameters (mu and nu)
 //' @export
 // [[Rcpp::export]]
 Rcpp::List beta_posterior_draw_cpp(const Rcpp::NumericVector& priorParams,
                                    double maxT_val, // Renamed maxT
                                    const Rcpp::NumericVector& mhStepSize_val, // Renamed mhStepSize
                                    const Rcpp::NumericMatrix& x_data, // Renamed x
                                    int n_draws = 1, int mhDrawsNum = 250) { // Renamed n and mhDraws to mhDrawsNum
   arma::mat x_arma = Rcpp::as<arma::mat>(x_data);
   return dp::BetaMixingDistribution::posteriorDrawStatic(priorParams, maxT_val, mhStepSize_val, x_arma, n_draws, mhDrawsNum);
 }

//' @title Beta distribution prior density (C++)
 //' @description C++ implementation for calculating the prior density for Beta parameters
 //' @param mu The mean parameter
 //' @param nu The precision parameter
 //' @param priorParams A numeric vector of prior parameters
 //' @param maxT The upper bound of the Beta distribution
 //' @return The prior density value
 //' @export
 // [[Rcpp::export]]
 double beta_prior_density_cpp(double mu_val, double nu_val, // Renamed params
                               const Rcpp::NumericVector& priorParams,
                               double maxT_val) { // Renamed maxT
   Rcpp::NumericVector mu_arr(1);
   Rcpp::NumericVector nu_arr(1);
   mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
   nu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
   mu_arr[0] = mu_val;
   nu_arr[0] = nu_val;
   Rcpp::List theta = Rcpp::List::create(mu_arr, nu_arr);

   dp::BetaMixingDistribution md(priorParams);
   md.maxT = maxT_val;
   return md.priorDensity(theta);
 }


//' @title Non-conjugate Beta Cluster Parameter Update (C++)
 //' @description C++ implementation for updating cluster parameters for a non-conjugate Beta DP model
 //' @param dpObj A list representing the Dirichlet Process object
 //' @return A list containing the updated cluster parameters
 //' @export
 // [[Rcpp::export]]
 Rcpp::List nonconjugate_beta_cluster_parameter_update_cpp(Rcpp::List dpObj) {
   // Extract necessary components from dpObj
   Rcpp::List mixingDistributionList = dpObj["mixingDistribution"];
   Rcpp::NumericVector priorParams = mixingDistributionList["priorParameters"];
   double maxT = Rcpp::as<double>(mixingDistributionList["maxT"]);
   Rcpp::NumericVector mhStepSize = Rcpp::as<Rcpp::NumericVector>(mixingDistributionList["mhStepSize"]);
   int mhDraws_val = Rcpp::as<int>(dpObj["mhDraws"]); // Get mhDraws from the dpObj and use it

   arma::mat data_mat = Rcpp::as<arma::mat>(dpObj["data"]); // Renamed data
   arma::uvec clusterLabels_vec = Rcpp::as<arma::uvec>(dpObj["clusterLabels"]); // Renamed clusterLabels
   int numberClusters_val = Rcpp::as<int>(dpObj["numberClusters"]); // Renamed numberClusters
   Rcpp::List clusterParameters_list = dpObj["clusterParameters"]; // Renamed clusterParameters


   Rcpp::NumericVector current_mus = Rcpp::as<Rcpp::NumericVector>(clusterParameters_list[0]);
   Rcpp::NumericVector current_nus = Rcpp::as<Rcpp::NumericVector>(clusterParameters_list[1]);

   dp::BetaMixingDistribution md(priorParams);
   md.maxT = maxT;
   md.mhStepSize = mhStepSize;

   for (int k = 0; k < numberClusters_val; ++k) {
     arma::uvec indices = arma::find(clusterLabels_vec == k);
     if (indices.n_elem > 0) {
       arma::mat cluster_data = data_mat.rows(indices);

       Rcpp::NumericVector mu_start_vec(1);
       Rcpp::NumericVector nu_start_vec(1);
       mu_start_vec.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
       nu_start_vec.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
       mu_start_vec[0] = current_mus[k];
       nu_start_vec[0] = current_nus[k];
       Rcpp::List start_pos = Rcpp::List::create(mu_start_vec, nu_start_vec);

       Rcpp::List mh_output = md.metropolisHastings(cluster_data, start_pos, mhDraws_val); // Use mhDraws_val
       Rcpp::List posterior_samples = mh_output["parameter_samples"];

       Rcpp::NumericVector mu_samples = posterior_samples[0];
       Rcpp::NumericVector nu_samples = posterior_samples[1];

       if(mu_samples.size() > 0) { // Check if MH produced samples
         current_mus[k] = mu_samples[mu_samples.size() - 1];  // Last sample
         current_nus[k] = nu_samples[nu_samples.size() - 1];
       }
     }
   }

   return Rcpp::List::create(
     Rcpp::Named("mu") = current_mus,
     Rcpp::Named("nu") = current_nus
   );
 }

//' @title Beta Metropolis-Hastings Sampler (C++)
 //' @description C++ implementation of a Metropolis-Hastings sampler for Beta distribution parameters
 //' @param x A numeric matrix of data points
 //' @param startMu Initial value for mu
 //' @param startNu Initial value for nu
 //' @param priorParams A numeric vector of prior parameters
 //' @param maxT The upper bound of the Beta distribution
 //' @param mhStepSize A numeric vector of step sizes for the MH algorithm
 //' @param noDraws The number of MH iterations
 //' @return A list containing the sampled parameters and acceptance ratio
 //' @export
 // [[Rcpp::export]]
 Rcpp::List beta_metropolis_hastings_cpp(const Rcpp::NumericMatrix& x_data, // Renamed x
                                         double startMu_val, double startNu_val, // Renamed params
                                         const Rcpp::NumericVector& priorParams,
                                         double maxT_val, // Renamed maxT
                                         const Rcpp::NumericVector& mhStepSize_val, // Renamed mhStepSize
                                         int noDraws_val = 100) { // Renamed noDraws
   // Create distribution object
   dp::BetaMixingDistribution md(priorParams);
   md.maxT = maxT_val;
   md.mhStepSize = mhStepSize_val;

   // Create start position
   Rcpp::NumericVector mu_start(1);
   Rcpp::NumericVector nu_start(1);
   mu_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
   nu_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
   mu_start[0] = startMu_val;
   nu_start[0] = startNu_val;

   Rcpp::List start_pos = Rcpp::List::create(mu_start, nu_start);

   // Convert data
   arma::mat x_arma = Rcpp::as<arma::mat>(x_data);

   // Run MH sampler directly as it's public
   return md.metropolisHastings(x_arma, start_pos, noDraws_val);
 }
