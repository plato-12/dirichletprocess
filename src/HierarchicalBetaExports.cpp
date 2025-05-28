// src/HierarchicalBetaExports.cpp
#include "../inst/include/HierarchicalDP.h"
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"

//' @title Create a Hierarchical Beta DP from R object (C++)
 //' @description C++ implementation for creating a Hierarchical Beta DP from an R object.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return An updated list with C++ object reference.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_create_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Store pointer as external pointer
   Rcpp::XPtr<dp::HierarchicalBetaDP> hdp_ptr(hdp, true);
   dpList.attr("cpp_ptr") = hdp_ptr;

   return dpList;
 }

//' @title Fit Hierarchical Beta DP (C++)
 //' @description C++ implementation for fitting a Hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @param iterations Number of iterations.
 //' @param updatePrior Whether to update prior parameters.
 //' @param progressBar Whether to show progress bar.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_fit_cpp(Rcpp::List dpList, int iterations,
                                      bool updatePrior = false,
                                      bool progressBar = true) {
   // Create C++ object from R
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Fit the model
   hdp->fit(iterations, updatePrior, progressBar);

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Add iteration info
   result["iterations"] = iterations;

   // Clean up
   delete hdp;

   return result;
 }

//' @title Update cluster components for Hierarchical Beta DP (C++)
 //' @description C++ implementation of cluster component update for hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_cluster_component_update_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Perform update
   hdp->clusterComponentUpdate();

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Clean up
   delete hdp;

   return result;
 }

//' @title Update cluster parameters for Hierarchical Beta DP (C++)
 //' @description C++ implementation of cluster parameter update for hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_cluster_parameter_update_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Perform update
   hdp->clusterParameterUpdate();

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Clean up
   delete hdp;

   return result;
 }

//' @title Update alpha for Hierarchical Beta DP (C++)
 //' @description C++ implementation of alpha update for hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_update_alpha_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Perform update
   hdp->updateAlpha();

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Clean up
   delete hdp;

   return result;
 }

//' @title Update global parameters for Hierarchical Beta DP (C++)
 //' @description C++ implementation of global parameter update for hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_global_parameter_update_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Perform update
   hdp->globalParameterUpdate();

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Clean up
   delete hdp;

   return result;
 }

//' @title Update G0 for Hierarchical Beta DP (C++)
 //' @description C++ implementation of G0 update for hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_update_g0_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Perform update
   hdp->updateG0();

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Clean up
   delete hdp;

   return result;
 }

//' @title Update gamma for Hierarchical Beta DP (C++)
 //' @description C++ implementation of gamma update for hierarchical Beta DP.
 //' @param dpList An R list representing the hierarchical DP object.
 //' @return Updated hierarchical DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_update_gamma_cpp(Rcpp::List dpList) {
   dp::HierarchicalBetaDP* hdp = dp::HierarchicalBetaDP::fromR(dpList);

   // Perform update
   hdp->updateGamma();

   // Convert back to R
   Rcpp::List result = hdp->toR();

   // Clean up
   delete hdp;

   return result;
 }

//' @title Create Hierarchical Beta mixing distributions (C++)
 //' @description C++ implementation for creating hierarchical Beta mixing distributions.
 //' @param n Number of datasets.
 //' @param priorParameters Prior parameters for the Beta distribution.
 //' @param hyperPriorParameters Hyper prior parameters.
 //' @param alphaPrior Alpha prior parameters.
 //' @param maxT Maximum value for Beta distribution.
 //' @param gammaPrior Gamma prior parameters.
 //' @param mhStepSize Metropolis-Hastings step size.
 //' @param num_sticks Number of stick breaking values.
 //' @return List of mixing distributions.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List hierarchical_beta_mixing_create_cpp(
     int n,
     Rcpp::NumericVector priorParameters,
     Rcpp::NumericVector hyperPriorParameters,
     Rcpp::NumericVector alphaPrior,
     double maxT,
     Rcpp::NumericVector gammaPrior,
     Rcpp::NumericVector mhStepSize,
     int num_sticks) {

   // Create base Beta mixing distribution
   dp::BetaMixingDistribution baseMD(priorParameters);
   baseMD.maxT = maxT;
   baseMD.mhStepSize = mhStepSize;
   baseMD.hyperPriorParameters = hyperPriorParameters;

   // Draw gamma
   double gamma = R::rgamma(gammaPrior[0], 1.0 / gammaPrior[1]);

   // Draw global parameters
   Rcpp::List theta_k = baseMD.priorDraw(num_sticks);

   // Stick breaking weights
   arma::vec beta_k(num_sticks);
   double remaining = 1.0;
   for (int i = 0; i < num_sticks - 1; i++) {
     double beta = R::rbeta(1.0, gamma);
     beta_k[i] = beta * remaining;
     remaining *= (1.0 - beta);
   }
   beta_k[num_sticks - 1] = remaining;

   // Create individual mixing distributions
   Rcpp::List mdobj_list(n);

   for (int i = 0; i < n; i++) {
     Rcpp::List mdobj;

     // Copy base distribution properties
     mdobj["distribution"] = "beta";
     mdobj["priorParameters"] = priorParameters;
     mdobj["conjugate"] = false;
     mdobj["mhStepSize"] = mhStepSize;
     mdobj["hyperPriorParameters"] = hyperPriorParameters;
     mdobj["maxT"] = maxT;

     // Hierarchical properties
     mdobj["theta_k"] = theta_k;
     mdobj["beta_k"] = Rcpp::wrap(beta_k);
     mdobj["gamma"] = gamma;

     // Individual alpha
     double alpha = R::rgamma(alphaPrior[0], 1.0 / alphaPrior[1]);
     mdobj["alpha"] = alpha;

     // Draw pi_k using stick breaking
     arma::vec pi_k(num_sticks);
     arma::vec beta_cumsum = arma::cumsum(beta_k);

     for (int j = 0; j < num_sticks; j++) {
       double shape2 = 1.0 - beta_cumsum[j];
       if (shape2 < 0) shape2 = 0;

       double pi_prime = R::rbeta(alpha * beta_k[j], alpha * shape2);

       // Compute stick breaking weight
       double prod = 1.0;
       for (int k = 0; k < j; k++) {
         prod *= (1.0 - pi_prime);
       }
       pi_k[j] = pi_prime * prod;
     }

     mdobj["pi_k"] = Rcpp::wrap(pi_k);

     // Set class
     mdobj.attr("class") = Rcpp::CharacterVector::create("hierarchical", "beta", "nonconjugate");

     mdobj_list[i] = mdobj;
   }

   return mdobj_list;
 }
