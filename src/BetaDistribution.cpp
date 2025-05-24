// src/BetaDistribution.cpp
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <cmath>

namespace dp {

// BetaMixingDistribution implementation
BetaMixingDistribution::BetaMixingDistribution(const Rcpp::NumericVector& priorParams) : maxT(1.0) {
  distribution = "beta";
  conjugate = false;
  priorParameters = priorParams;

  // Set default values for mhStepSize and hyperPriorParameters if not provided
  mhStepSize = Rcpp::NumericVector::create(1.0, 1.0);
  hyperPriorParameters = Rcpp::NumericVector::create(1.0, 0.125);
}

BetaMixingDistribution::~BetaMixingDistribution() {
  // Destructor
}

Rcpp::NumericVector BetaMixingDistribution::likelihood(const arma::vec& x_data, const Rcpp::List& theta) const { // Renamed x to x_data
  // Extract parameters
  Rcpp::NumericVector mu_array = theta[0];
  Rcpp::NumericVector nu_array = theta[1];

  int n_data = x_data.n_elem;
  Rcpp::NumericVector result(n_data);

  // For each cluster parameter
  double mu = mu_array[0];
  double tau = nu_array[0];

  // Transform to Beta parameters
  double a = (mu * tau) / maxT;
  double b = (1.0 - mu/maxT) * tau;

  // Calculate likelihood for each data point
  for (int i = 0; i < n_data; i++) {
    if (x_data[i] >= 0 && x_data[i] <= maxT) {
      result[i] = (1.0/maxT) * R::dbeta(x_data[i]/maxT, a, b, false);
    } else {
      result[i] = 0.0;
    }
  }

  return result;
}

Rcpp::List BetaMixingDistribution::priorDraw(int n_draws) const { // Renamed n to n_draws

  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  Rcpp::Rcout << "C++ priorParams: " << priorParams[0] << " " << priorParams[1] << std::endl;
  Rcpp::Rcout << "C++ maxT: " << maxT << std::endl;

  Rcpp::NumericVector mu(n_draws);
  Rcpp::NumericVector nu(n_draws);

  for (int i = 0; i < n_draws; i++) {
    // mu ~ Uniform(0, maxT)
    mu[i] = R::runif(0.0, maxT);

    // nu ~ InverseGamma(priorParams[0], priorParams[1])
    // nu = 1/gamma where gamma ~ Gamma(priorParams[0], 1/priorParams[1])
    double gamma_draw = R::rgamma(priorParams[0], 1.0/priorParams[1]);
    nu[i] = 1.0 / gamma_draw;
  }

  // Convert to 3D arrays with dimension (1,1,n_draws)
  Rcpp::NumericVector mu_arr(n_draws);
  Rcpp::NumericVector nu_arr(n_draws);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n_draws);
  nu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n_draws);

  for (int i = 0; i < n_draws; i++) {
    mu_arr[i] = mu[i];
    nu_arr[i] = nu[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("nu") = nu_arr
  );
}

double BetaMixingDistribution::priorDensity(const Rcpp::List& theta) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  Rcpp::NumericVector mu_array = theta[0];
  Rcpp::NumericVector nu_array = theta[1];

  double mu = mu_array[0];
  double nu = nu_array[0];

  // mu ~ Uniform(0, maxT)
  double muDensity = (mu >= 0 && mu <= maxT) ? 1.0/maxT : 0.0;

  // nu ~ InverseGamma(priorParams[0], priorParams[1])
  // Using the relationship: if X ~ IG(a,b), then pdf(x) = (b^a/Gamma(a)) * x^(-a-1) * exp(-b/x)
  double nuDensity = 0.0;
  if (nu > 0) {
    double shape = priorParams[0];
    double scale = priorParams[1];
    nuDensity = R::dgamma(1.0/nu, shape, 1.0/scale, false) / (nu * nu);
  }

  return muDensity * nuDensity;
}

Rcpp::List BetaMixingDistribution::mhParameterProposal(const Rcpp::List& oldParams) const {
  Rcpp::NumericVector mhStep = Rcpp::as<Rcpp::NumericVector>(mhStepSize);

  Rcpp::NumericVector old_mu = oldParams[0];
  Rcpp::NumericVector old_nu = oldParams[1];

  double mu = old_mu[0];
  double nu = old_nu[0];

  // Propose new mu
  double new_mu = mu + mhStep[0] * R::rnorm(0.0, 2.4);

  // Reflect at boundaries
  if (new_mu < 0 || new_mu > maxT) {
    new_mu = mu;  // Reject proposals outside bounds
  }

  // Propose new nu (must be positive)
  double new_nu = std::abs(nu + mhStep[1] * R::rnorm(0.0, 2.4));

  // Create return arrays
  Rcpp::NumericVector mu_arr(1);
  Rcpp::NumericVector nu_arr(1);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  nu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

  mu_arr[0] = new_mu;
  nu_arr[0] = new_nu;

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("nu") = nu_arr
  );
}

Rcpp::List BetaMixingDistribution::metropolisHastings(const arma::mat& x_data, // Renamed x to x_data
                                                      const Rcpp::List& startPos,
                                                      int noDraws) const {
  // Initialize parameter storage
  Rcpp::NumericVector mu_samples(noDraws);
  Rcpp::NumericVector nu_samples(noDraws);

  // Set initial values
  Rcpp::NumericVector start_mu = startPos[0];
  Rcpp::NumericVector start_nu = startPos[1];

  Rcpp::List current_params = Rcpp::clone(startPos);

  // Calculate initial log-likelihood and log-prior
  double current_log_lik = 0.0;
  Rcpp::NumericVector lik_vals = likelihood(arma::vectorise(x_data), current_params);
  for (int i = 0; i < lik_vals.size(); i++) {
    if (lik_vals[i] > 0) {
      current_log_lik += std::log(lik_vals[i]);
    } else {
      current_log_lik += -1e10;  // Large negative value for numerical stability
    }
  }
  double current_log_prior = std::log(priorDensity(current_params));

  int accept_count = 0;

  // Store first sample
  mu_samples[0] = start_mu[0];
  nu_samples[0] = start_nu[0];

  // Main MCMC loop
  for (int iter = 1; iter < noDraws; iter++) {
    // Propose new parameters
    Rcpp::List proposed_params = mhParameterProposal(current_params);

    // Calculate proposed log-likelihood
    double proposed_log_lik = 0.0;
    Rcpp::NumericVector prop_lik_vals = likelihood(arma::vectorise(x_data), proposed_params);
    for (int i = 0; i < prop_lik_vals.size(); i++) {
      if (prop_lik_vals[i] > 0) {
        proposed_log_lik += std::log(prop_lik_vals[i]);
      } else {
        proposed_log_lik += -1e10;
      }
    }

    // Calculate proposed log-prior
    double proposed_log_prior = std::log(priorDensity(proposed_params));

    // Calculate acceptance ratio
    double log_ratio = (proposed_log_lik + proposed_log_prior) -
      (current_log_lik + current_log_prior);
    double accept_prob = std::min(1.0, std::exp(log_ratio));

    // Accept or reject
    if (R::runif(0, 1) < accept_prob) {
      current_params = proposed_params;
      current_log_lik = proposed_log_lik;
      current_log_prior = proposed_log_prior;
      accept_count++;
    }

    // Store current values
    Rcpp::NumericVector curr_mu = current_params[0];
    Rcpp::NumericVector curr_nu = current_params[1];
    mu_samples[iter] = curr_mu[0];
    nu_samples[iter] = curr_nu[0];
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_samples,
    Rcpp::Named("nu") = nu_samples
  );
}

Rcpp::List BetaMixingDistribution::posteriorDraw(const arma::mat& x_data, int n_draws) const {
  // Start from a prior draw
  Rcpp::List start_pos = priorDraw(1);

  // Run Metropolis-Hastings
  Rcpp::List mh_result = metropolisHastings(x_data, start_pos, n_draws);

  return mh_result;
}

void BetaMixingDistribution::updatePriorParameters(const Rcpp::List& clusterParameters, int n_clusters) {
  Rcpp::NumericVector hyperPrior = Rcpp::as<Rcpp::NumericVector>(hyperPriorParameters);
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  Rcpp::NumericVector nu_values = clusterParameters[1];
  int num_actual_clusters = nu_values.size();

  double posteriorShape = hyperPrior[0] + priorParams[0] * num_actual_clusters;
  double posteriorRate = hyperPrior[1];
  if(num_actual_clusters > 0) {
    posteriorRate += arma::sum(1.0 / Rcpp::as<arma::vec>(nu_values));
  }

  double newGamma = R::rgamma(posteriorShape, 1.0/posteriorRate);

  Rcpp::NumericVector newPriorParams = Rcpp::NumericVector::create(priorParams[0], newGamma);
  priorParameters = newPriorParams;
}

// Static methods for direct testing
Rcpp::List BetaMixingDistribution::priorDrawStatic(const Rcpp::NumericVector& priorParams,
                                                   double maxT_val, int n_draws) {
  BetaMixingDistribution md(priorParams);
  md.maxT = maxT_val;
  return md.priorDraw(n_draws);
}

Rcpp::List BetaMixingDistribution::posteriorDrawStatic(const Rcpp::NumericVector& priorParams,
                                                       double maxT_val,
                                                       const Rcpp::NumericVector& mhStepSize_val,
                                                       const arma::mat& x_data,
                                                       int n_draws, int mhDrawsVal) {
  BetaMixingDistribution md(priorParams);
  md.maxT = maxT_val;
  md.mhStepSize = mhStepSize_val;

  Rcpp::List start_pos = md.priorDraw(1);
  Rcpp::List mh_result = md.metropolisHastings(x_data, start_pos, mhDrawsVal);

  Rcpp::NumericVector mu_all = mh_result["mu"];
  Rcpp::NumericVector nu_all = mh_result["nu"];

  Rcpp::NumericVector mu_final(n_draws);
  Rcpp::NumericVector nu_final(n_draws);
  mu_final.attr("dim") = Rcpp::IntegerVector::create(1, 1, n_draws);
  nu_final.attr("dim") = Rcpp::IntegerVector::create(1, 1, n_draws);

  if (mhDrawsVal >= n_draws) {
    int start_idx = mhDrawsVal - n_draws;
    for (int i = 0; i < n_draws; i++) {
      mu_final[i] = mu_all[start_idx + i];
      nu_final[i] = nu_all[start_idx + i];
    }
  } else {
    for (int i = 0; i < mhDrawsVal; i++) {
      mu_final[i] = mu_all[i];
      nu_final[i] = nu_all[i];
    }
    for (int i = mhDrawsVal; i < n_draws; i++) {
      mu_final[i] = mu_all[mhDrawsVal > 0 ? mhDrawsVal -1 : 0];
      nu_final[i] = nu_all[mhDrawsVal > 0 ? mhDrawsVal -1 : 0];
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_final,
    Rcpp::Named("nu") = nu_final
  );
}

Rcpp::NumericVector BetaMixingDistribution::likelihoodStatic(const arma::vec& x_data,
                                                             double mu_val, double nu_val,
                                                             double maxT_val) {
  Rcpp::NumericVector mu_arr(1);
  Rcpp::NumericVector nu_arr(1);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  nu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  mu_arr[0] = mu_val;
  nu_arr[0] = nu_val;

  Rcpp::List theta = Rcpp::List::create(mu_arr, nu_arr);

  BetaMixingDistribution md(Rcpp::NumericVector::create(2, 8));
  md.maxT = maxT_val;

  return md.likelihood(x_data, theta);
}

// NonConjugateBetaDP implementation
// Constructor and Destructor implementations should be in BetaDP.cpp
// NonConjugateBetaDP::NonConjugateBetaDP() : mixingDistribution(nullptr), numberClusters(0), m(3) {
// }
// NonConjugateBetaDP::~NonConjugateBetaDP() {
//   if (mixingDistribution) {
//     delete mixingDistribution;
//   }
// }

// Implementations for clusterComponentUpdate and clusterLabelChange
// should be *only* in BetaDP.cpp to avoid multiple definition errors.
// void NonConjugateBetaDP::clusterComponentUpdate() {
//   Rcpp::stop("NonConjugateBetaDP::clusterComponentUpdate - This definition should be removed from BetaDistribution.cpp");
// }

// Rcpp::List NonConjugateBetaDP::clusterLabelChange(int i, int newLabel, int currentLabel,
//                                                   const Rcpp::List& aux) {
//   Rcpp::stop("NonConjugateBetaDP::clusterLabelChange - This definition should be removed from BetaDistribution.cpp");
// }


// Corrected clusterParameterUpdate for NonConjugateBetaDP
void NonConjugateBetaDP::clusterParameterUpdate() {
  for (int k = 0; k < numberClusters; k++) {
    arma::uvec clusterIndices = arma::find(clusterLabels == k);

    if (clusterIndices.n_elem > 0) {
      arma::mat clusterData = data.rows(clusterIndices);

      Rcpp::NumericVector mu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector nu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      Rcpp::NumericVector mu_start(1);
      Rcpp::NumericVector nu_start(1);
      mu_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      nu_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      mu_start[0] = mu_vec[k];
      nu_start[0] = nu_vec[k];

      Rcpp::List start_pos = Rcpp::List::create(mu_start, nu_start);

      int n_mh_draws = 250; // Default MH draws.
      // This should ideally be a member of NonConjugateBetaDP,
      // set from the R dpObj's mhDraws parameter.
      // The R-exported function in BetaExports.cpp handles passing
      // the dpObj[["mhDraws"]] to BetaMixingDistribution::metropolisHastings.

      Rcpp::List mh_result = mixingDistribution->metropolisHastings(clusterData, start_pos, n_mh_draws);

      Rcpp::NumericVector mu_samples = mh_result["mu"];
      Rcpp::NumericVector nu_samples = mh_result["nu"];

      if(mu_samples.size() > 0) {
        mu_vec[k] = mu_samples[mu_samples.size() - 1];
        nu_vec[k] = nu_samples[nu_samples.size() - 1];
      }

      clusterParameters[0] = mu_vec;
      clusterParameters[1] = nu_vec;
    }
  }
}

void NonConjugateBetaDP::updateAlpha() {
  double x_draw_val = R::rbeta(alpha + 1.0, n);

  Rcpp::NumericVector alphaPriors = Rcpp::as<Rcpp::NumericVector>(alphaPriorParameters);

  double pi1 = alphaPriors[0] + numberClusters - 1.0;
  double pi2 = n * (alphaPriors[1] - log(x_draw_val));
  double pi_ratio = pi1 / (pi1 + pi2);

  double postShape, postRate;
  if (R::runif(0, 1) < pi_ratio) {
    postShape = alphaPriors[0] + numberClusters;
  } else {
    postShape = alphaPriors[0] + numberClusters - 1.0;
  }
  postRate = alphaPriors[1] - log(x_draw_val);

  alpha = R::rgamma(postShape, 1.0/postRate);
}

} // namespace dp
