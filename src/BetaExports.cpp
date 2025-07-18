// src/BetaExports.cpp

#include <RcppArmadillo.h>
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/DirichletProcess.h"
#include <iostream>

// Helper function to get a single data point (row)
arma::rowvec get_row(const arma::mat& m, int i) {
  return m.row(i);
}

// [[Rcpp::export]]
Rcpp::List beta_prior_draw_cpp(const Rcpp::NumericVector& priorParams, double maxT, int n) {
  return dp::BetaMixingDistribution::priorDrawStatic(priorParams, maxT, n);
}

// [[Rcpp::export]]
Rcpp::NumericVector beta_likelihood_cpp(const Rcpp::NumericVector& x, double mu, double nu, double maxT) {
  arma::vec x_arma = Rcpp::as<arma::vec>(x);
  return dp::BetaMixingDistribution::likelihoodStatic(x_arma, mu, nu, maxT);
}

// [[Rcpp::export]]
double beta_prior_density_cpp(double mu, double nu, const Rcpp::NumericVector& priorParams, double maxT) {
  dp::BetaMixingDistribution md(priorParams);
  md.maxT = maxT;

  Rcpp::NumericVector mu_arr(1, mu);
  Rcpp::NumericVector nu_arr(1, nu);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  nu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

  Rcpp::List theta = Rcpp::List::create(mu_arr, nu_arr);
  return md.priorDensity(theta);
}

// [[Rcpp::export]]
Rcpp::List beta_metropolis_hastings_cpp(const Rcpp::NumericMatrix& x, double startMu, double startNu,
                                        const Rcpp::NumericVector& priorParams, double maxT,
                                        const Rcpp::NumericVector& mhStep, int noDraws) {
  dp::BetaMixingDistribution md(priorParams);
  md.maxT = maxT;
  md.mhStepSize = mhStep;

  Rcpp::NumericVector mu_start(1, startMu);
  Rcpp::NumericVector nu_start(1, startNu);
  mu_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  nu_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

  Rcpp::List startPos = Rcpp::List::create(mu_start, nu_start);
  arma::mat x_arma = Rcpp::as<arma::mat>(x);
  return md.metropolisHastings(x_arma, startPos, noDraws);
}

// [[Rcpp::export]]
Rcpp::List beta_posterior_draw_cpp(const Rcpp::NumericVector& priorParams, double maxT_val,
                                   const Rcpp::NumericVector& mhStepSize_val, const Rcpp::NumericMatrix& x_data,
                                   int n_draws, int mhDrawsVal) {
  arma::mat x_arma = Rcpp::as<arma::mat>(x_data);
  return dp::BetaMixingDistribution::posteriorDrawStatic(priorParams, maxT_val, mhStepSize_val, x_arma, n_draws, mhDrawsVal);
}

// [[Rcpp::export]]
SEXP nonconjugate_beta_cluster_parameter_update_cpp(Rcpp::List dp_list) {
  Rcpp::warning("C++ function 'nonconjugate_beta_cluster_parameter_update_cpp' is a STUB and not implemented. R fallback should be used.");
  return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::List nonconjugate_beta_cluster_component_update_cpp(Rcpp::List dp_list) {
  try {
    // Extract necessary components from R list
    arma::mat data = Rcpp::as<arma::mat>(dp_list["data"]);
    arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dp_list["clusterLabels"]);
    arma::uvec pointsPerCluster = Rcpp::as<arma::uvec>(dp_list["pointsPerCluster"]);
    int numberClusters = dp_list["numberClusters"];
    double alpha = dp_list["alpha"];
    Rcpp::List mixingDistribution = dp_list["mixingDistribution"];
    Rcpp::List clusterParameters = dp_list["clusterParameters"];
    int m = dp_list.containsElementNamed("m") ? Rcpp::as<int>(dp_list["m"]) : 3;
    
    // Convert R's 1-based indexing to C++'s 0-based indexing
    clusterLabels = clusterLabels - 1;
    
    // Extract Beta-specific parameters
    Rcpp::NumericVector priorParams = mixingDistribution["priorParameters"];
    double maxT = mixingDistribution.containsElementNamed("maxT") ? 
                  Rcpp::as<double>(mixingDistribution["maxT"]) : 1.0;
    
    // Create C++ DP object
    std::unique_ptr<dp::NonConjugateBetaDP> dp_cpp(new dp::NonConjugateBetaDP());
    
    // Initialize the C++ object
    dp_cpp->data = data;
    dp_cpp->n = data.n_rows;
    dp_cpp->alpha = alpha;
    dp_cpp->clusterLabels = clusterLabels;
    dp_cpp->pointsPerCluster = pointsPerCluster;
    dp_cpp->numberClusters = numberClusters;
    dp_cpp->clusterParameters = clusterParameters;
    dp_cpp->m = m;
    
    // Create and initialize mixing distribution
    dp_cpp->mixingDistribution = std::unique_ptr<dp::BetaMixingDistribution>(
      new dp::BetaMixingDistribution(priorParams));
    dp_cpp->mixingDistribution->maxT = maxT;
    
    // Copy mhStepSize if available
    if (mixingDistribution.containsElementNamed("mhStepSize")) {
      dp_cpp->mixingDistribution->mhStepSize = mixingDistribution["mhStepSize"];
    }
    
    // Copy hyperPriorParameters if available
    if (mixingDistribution.containsElementNamed("hyperPriorParameters")) {
      dp_cpp->mixingDistribution->hyperPriorParameters = mixingDistribution["hyperPriorParameters"];
    }
    
    // Perform cluster component update
    dp_cpp->clusterComponentUpdate();
    
    // Convert C++'s 0-based indexing back to R's 1-based indexing
    arma::uvec clusterLabels_R = dp_cpp->clusterLabels + 1;
    
    // Return updated R list (keep original structure intact)
    Rcpp::List result = Rcpp::clone(dp_list);
    result["clusterLabels"] = clusterLabels_R;
    result["pointsPerCluster"] = dp_cpp->pointsPerCluster;
    result["numberClusters"] = dp_cpp->numberClusters;
    result["clusterParameters"] = dp_cpp->clusterParameters;
    
    return result;
    
  } catch (const std::exception& e) {
    Rcpp::stop("Error in nonconjugate_beta_cluster_component_update_cpp: " + std::string(e.what()));
  }
}
