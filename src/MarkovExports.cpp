// src/MarkovExports.cpp
#include "../inst/include/MarkovDP.h"
#include "../inst/include/RcppConversions.h"

//' @title Create a Markov DP from R object (C++)
 //' @description C++ implementation for creating a Markov DP from an R object.
 //' @param dpObj An R list representing the Markov DP object.
 //' @return An updated list with C++ object reference.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List markov_dp_create_cpp(Rcpp::List dpObj) {
   dp::MarkovDP* mdp = dp::MarkovDP::fromR(dpObj);

   // Store pointer as external pointer
   Rcpp::XPtr<dp::MarkovDP> mdp_ptr(mdp, true);
   dpObj.attr("cpp_ptr") = mdp_ptr;

   return dpObj;
 }

//' @title Fit Markov DP (C++)
 //' @description C++ implementation for fitting a Markov DP (HMM).
 //' @param dpObj An R list representing the Markov DP object.
 //' @param iterations Number of iterations.
 //' @param updatePrior Whether to update prior parameters.
 //' @param progressBar Whether to show progress bar.
 //' @return Updated Markov DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List markov_dp_fit_cpp(Rcpp::List dpObj, int iterations,
                              bool updatePrior = false,
                              bool progressBar = true) {
   // Create C++ object from R
   dp::MarkovDP* mdp = dp::MarkovDP::fromR(dpObj);

   // Fit the model
   mdp->fit(iterations, updatePrior, progressBar);

   // Convert back to R
   Rcpp::List result = mdp->toR();

   // Clean up
   delete mdp;

   return result;
 }

//' @title Update states for Markov DP (C++)
 //' @description C++ implementation of state update for Markov DP.
 //' @param dpObj An R list representing the Markov DP object.
 //' @return Updated Markov DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List markov_dp_update_states_cpp(Rcpp::List dpObj) {
   dp::MarkovDP* mdp = dp::MarkovDP::fromR(dpObj);

   // Perform update
   mdp->updateStates();

   // Convert back to R
   Rcpp::List result = mdp->toR();

   // Clean up
   delete mdp;

   return result;
 }

//' @title Update alpha and beta for Markov DP (C++)
 //' @description C++ implementation of alpha/beta update for Markov DP.
 //' @param dpObj An R list representing the Markov DP object.
 //' @return Updated Markov DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List markov_dp_update_alpha_beta_cpp(Rcpp::List dpObj) {
   dp::MarkovDP* mdp = dp::MarkovDP::fromR(dpObj);

   // Perform update
   mdp->updateAlphaBeta();

   // Convert back to R
   Rcpp::List result = mdp->toR();

   // Clean up
   delete mdp;

   return result;
 }

//' @title Update parameters for Markov DP (C++)
 //' @description C++ implementation of parameter update for Markov DP.
 //' @param dpObj An R list representing the Markov DP object.
 //' @return Updated Markov DP object.
 //' @export
 // [[Rcpp::export]]
 Rcpp::List markov_dp_param_update_cpp(Rcpp::List dpObj) {
   dp::MarkovDP* mdp = dp::MarkovDP::fromR(dpObj);

   // Perform update
   mdp->paramUpdate();

   // Convert back to R
   Rcpp::List result = mdp->toR();

   // Clean up
   delete mdp;

   return result;
 }
