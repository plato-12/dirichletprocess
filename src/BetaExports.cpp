// src/BetaExports.cpp

#include <RcppArmadillo.h>
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/DirichletProcess.h"
#include <iostream> // For Rcpp::Rcout

// Helper function to get a single data point (row)
arma::rowvec get_row(const arma::mat& m, int i) {
  return m.row(i);
}

// [[Rcpp::export]]
Rcpp::List beta_prior_draw_cpp(const Rcpp::NumericVector& priorParams, double maxT, int n) {
  Rcpp::Rcout << "C++ priorParams: " << priorParams[0] << " " << priorParams[1] << std::endl;
  Rcpp::Rcout << "C++ maxT: " << maxT << std::endl;
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
Rcpp::List nonconjugate_beta_cluster_parameter_update_cpp(Rcpp::List dp_list) {
  Rcpp::warning("C++ function 'nonconjugate_beta_cluster_parameter_update_cpp' is not implemented.");
  return dp_list;
}

// [[Rcpp::export]]
Rcpp::List nonconjugate_beta_cluster_component_update_cpp(Rcpp::List dp_list) {

  arma::mat data = Rcpp::as<arma::mat>(dp_list["data"]);
  int n = data.n_rows;
  double alpha = Rcpp::as<double>(dp_list["alpha"]);
  arma::uvec clusterLabels = Rcpp::as<arma::uvec>(dp_list["clusterLabels"]);
  Rcpp::List clusterParameters = Rcpp::as<Rcpp::List>(dp_list["clusterParameters"]);
  int numberClusters = Rcpp::as<int>(dp_list["numberClusters"]);
  arma::uvec pointsPerCluster = Rcpp::as<arma::uvec>(dp_list["pointsPerCluster"]);
  int m = Rcpp::as<int>(dp_list["m"]);

  Rcpp::List md_list = Rcpp::as<Rcpp::List>(dp_list["mixingDistribution"]);
  Rcpp::NumericVector priorParameters_md = Rcpp::as<Rcpp::NumericVector>(md_list["priorParameters"]); // Renamed to avoid conflict
  dp::BetaMixingDistribution md(priorParameters_md);
  md.maxT = Rcpp::as<double>(md_list["maxT"]);
  md.mhStepSize = Rcpp::as<Rcpp::NumericVector>(md_list["mhStepSize"]);

  for (int i = 0; i < n; ++i) {

    int currentLabel = clusterLabels[i];
    pointsPerCluster[currentLabel]--;

    Rcpp::List auxParams;
    if (pointsPerCluster[currentLabel] == 0) {
      auxParams = md.priorDraw(m - 1); // This returns a List with "mu" and "nu"
      Rcpp::NumericVector current_mu_vec = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(clusterParameters)[0]);
      Rcpp::NumericVector current_nu_vec = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(clusterParameters)[1]);

      Rcpp::NumericVector aux_mu_draws = Rcpp::as<Rcpp::NumericVector>(auxParams["mu"]);
      Rcpp::NumericVector aux_nu_draws = Rcpp::as<Rcpp::NumericVector>(auxParams["nu"]);

      Rcpp::NumericVector combined_mu(m);
      Rcpp::NumericVector combined_nu(m);

      combined_mu[0] = current_mu_vec[currentLabel];
      combined_nu[0] = current_nu_vec[currentLabel];
      for(int k=0; k < m-1; ++k) {
        combined_mu[k+1] = aux_mu_draws[k]; // priorDraw returns 3D array, access as flat vector
        combined_nu[k+1] = aux_nu_draws[k];
      }
      combined_mu.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);
      combined_nu.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);
      auxParams = Rcpp::List::create(Rcpp::Named("mu") = combined_mu, Rcpp::Named("nu") = combined_nu);
    } else {
      auxParams = md.priorDraw(m); // This returns a List with "mu" and "nu"
    }

    Rcpp::NumericVector probs(numberClusters + m);
    arma::rowvec y_i = get_row(data, i);

    Rcpp::NumericVector current_cluster_mus = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(clusterParameters)[0]);
    Rcpp::NumericVector current_cluster_nus = Rcpp::as<Rcpp::NumericVector>(Rcpp::as<Rcpp::List>(clusterParameters)[1]);

    for (int j = 0; j < numberClusters; ++j) {
      if (pointsPerCluster[j] > 0) {
        Rcpp::List theta_j = Rcpp::List::create(
          Rcpp::Named("mu") = Rcpp::NumericVector::create(current_cluster_mus[j]),
          Rcpp::Named("nu") = Rcpp::NumericVector::create(current_cluster_nus[j])
        );
        Rcpp::as<Rcpp::NumericVector>(theta_j["mu"]).attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        Rcpp::as<Rcpp::NumericVector>(theta_j["nu"]).attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

        probs[j] = pointsPerCluster[j] * md.likelihood(y_i, theta_j)[0];
      } else {
        probs[j] = 0;
      }
    }

    Rcpp::NumericVector aux_mu_samples = Rcpp::as<Rcpp::NumericVector>(auxParams["mu"]);
    Rcpp::NumericVector aux_nu_samples = Rcpp::as<Rcpp::NumericVector>(auxParams["nu"]);
    for (int j = 0; j < m; ++j) {
      Rcpp::List theta_aux_j = Rcpp::List::create(
        Rcpp::Named("mu") = Rcpp::NumericVector::create(aux_mu_samples[j]), // Access as flat vector
        Rcpp::Named("nu") = Rcpp::NumericVector::create(aux_nu_samples[j])  // Access as flat vector
      );
      Rcpp::as<Rcpp::NumericVector>(theta_aux_j["mu"]).attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      Rcpp::as<Rcpp::NumericVector>(theta_aux_j["nu"]).attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

      probs[numberClusters + j] = (alpha / m) * md.likelihood(y_i, theta_aux_j)[0];
    }

    for(int k=0; k < probs.size(); ++k){
      if(R_IsNA(probs[k]) || !R_finite(probs[k])) probs[k] = 0;
    }
    if(Rcpp::sum(probs) == 0.0) {
      std::fill(probs.begin(), probs.end(), 1.0);
    }

    // CORRECTED: Using Rcpp::sample correctly
    Rcpp::IntegerVector newLabel_vec = Rcpp::sample(probs.size(), 1, true, probs);
    int newLabel_1based = newLabel_vec[0];

    clusterLabels[i] = newLabel_1based -1;
    pointsPerCluster[currentLabel]++;
  }

  dp_list["clusterLabels"] = clusterLabels;
  dp_list["pointsPerCluster"] = pointsPerCluster;

  return dp_list;
}
