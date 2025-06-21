// src/WeibullDistribution.cpp
#include "../inst/include/WeibullDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <cmath>
#include <limits>
#include <memory>

namespace dp {

// WeibullMixingDistribution implementation
WeibullMixingDistribution::WeibullMixingDistribution(const Rcpp::NumericVector& priorParams,
                                                     const Rcpp::NumericVector& mhStepSize,
                                                     const Rcpp::NumericVector& hyperPriorParams) {
  distribution = "weibull";
  conjugate = false;
  priorParameters = priorParams;
  this->mhStepSize = mhStepSize;

  if (hyperPriorParams.size() > 0) {
    hyperPriorParameters = hyperPriorParams;
  } else {
    hyperPriorParameters = Rcpp::NumericVector::create(6.0, 2.0, 1.0, 0.5);
  }
}

WeibullMixingDistribution::~WeibullMixingDistribution() {
  // Destructor
}

Rcpp::NumericVector WeibullMixingDistribution::likelihood(const arma::vec& x, const Rcpp::List& theta) const {
  Rcpp::NumericVector alpha_array = theta[0];
  Rcpp::NumericVector lambda_array = theta[1];

  int n_data = x.n_elem;
  Rcpp::NumericVector result(n_data);

  double alpha = alpha_array[0];
  double lambda = lambda_array[0];

  for (int i = 0; i < n_data; i++) {
    if (x[i] < 0) {
      result[i] = 0.0;  // Weibull not defined for negative values
    } else if (x[i] == 0) {
      // Special case for x = 0
      result[i] = (alpha == 1.0) ? lambda : 0.0;
    } else if (lambda > 0 && alpha > 0 && std::isfinite(lambda) && std::isfinite(alpha)) {
      // Use log-space calculation for stability
      double log_lik = std::log(alpha) - std::log(lambda) +
        (alpha - 1.0) * std::log(x[i]) -
        std::pow(x[i] / std::pow(lambda, 1.0/alpha), alpha);

      if (std::isfinite(log_lik) && log_lik > -700) {  // Prevent underflow
        result[i] = std::exp(log_lik);
      } else {
        result[i] = 1e-300;
      }
    } else {
      result[i] = 1e-300;
    }
  }

  return result;
}

Rcpp::List WeibullMixingDistribution::priorDraw(int n) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  Rcpp::NumericVector alpha_values(n);
  Rcpp::NumericVector lambda_values(n);

  // priorParams[0] = phi, priorParams[1] = alpha0, priorParams[2] = beta0
  for (int i = 0; i < n; i++) {
    alpha_values[i] = R::runif(0.0, priorParams[0]);
    // R's rgamma uses shape and scale, convert rate to scale: scale = 1/rate
    double gamma_draw = R::rgamma(priorParams[1], 1.0 / priorParams[2]);
    lambda_values[i] = 1.0 / gamma_draw;
  }

  // Convert to 3D arrays
  Rcpp::NumericVector alpha_arr(n);
  Rcpp::NumericVector lambda_arr(n);
  alpha_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);
  lambda_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    alpha_arr[i] = alpha_values[i];
    lambda_arr[i] = lambda_values[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha_arr,
    Rcpp::Named("lambda") = lambda_arr
  );
}

Rcpp::NumericVector WeibullMixingDistribution::priorDensity(const Rcpp::List& theta) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);
  Rcpp::NumericVector alpha_array = theta[0];

  double alpha = alpha_array[0];
  double density = 0.0;

  if (alpha > 0 && alpha < priorParams[0]) {
    density = 1.0 / priorParams[0];
  } else {
    density = 0.0;
  }

  return Rcpp::NumericVector::create(density);
}

Rcpp::List WeibullMixingDistribution::mhParameterProposal(const Rcpp::List& oldParams) const {
  Rcpp::NumericVector mhStep = Rcpp::as<Rcpp::NumericVector>(mhStepSize);
  Rcpp::NumericVector old_alpha = oldParams[0];
  Rcpp::NumericVector old_lambda = oldParams[1];

  double alpha_old = old_alpha[0];
  double new_alpha = std::abs(alpha_old + mhStep[0] * R::rnorm(0.0, 1.7));

  // Create return arrays
  Rcpp::NumericVector alpha_arr(1);
  Rcpp::NumericVector lambda_arr(1);
  alpha_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  lambda_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

  alpha_arr[0] = new_alpha;
  lambda_arr[0] = old_lambda[0]; // Lambda will be updated analytically later

  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha_arr,
    Rcpp::Named("lambda") = lambda_arr
  );
}

Rcpp::List WeibullMixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);
  Rcpp::NumericVector mhStep = Rcpp::as<Rcpp::NumericVector>(mhStepSize);

  int mhDraws = std::max(100, n * 10);  // Ensure enough samples
  const int MAX_ITER = 10000;  // Maximum iterations to prevent infinite loops

  // Initialize from prior
  Rcpp::List initial_draw = priorDraw(1);
  Rcpp::NumericVector alpha_init = initial_draw[0];
  Rcpp::NumericVector lambda_init = initial_draw[1];

  double alpha_current = alpha_init[0];
  double lambda_current = lambda_init[0];

  // Storage for samples
  std::vector<double> alpha_samples;
  std::vector<double> lambda_samples;
  alpha_samples.reserve(mhDraws);
  lambda_samples.reserve(mhDraws);

  // Create initial parameter list
  Rcpp::NumericVector alpha_vec(1), lambda_vec(1);
  alpha_vec.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  lambda_vec.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  alpha_vec[0] = alpha_current;
  lambda_vec[0] = lambda_current;
  Rcpp::List current_params = Rcpp::List::create(alpha_vec, lambda_vec);

  // Calculate initial likelihood and prior
  double current_log_lik = 0.0;
  Rcpp::NumericVector lik_vals = likelihood(arma::vectorise(x), current_params);
  for (int k = 0; k < lik_vals.size(); k++) {
    if (lik_vals[k] > 1e-300) {
      current_log_lik += std::log(lik_vals[k]);
    } else {
      current_log_lik = -std::numeric_limits<double>::infinity();
      break;
    }
  }
  double current_log_prior = std::log(priorDensity(current_params)[0]);

  int accept_count = 0;
  int iter_count = 0;

  for (int iter = 0; iter < mhDraws && iter_count < MAX_ITER; iter++, iter_count++) {
    // Propose new alpha
    Rcpp::List proposed_params = mhParameterProposal(current_params);
    Rcpp::NumericVector prop_alpha = proposed_params[0];
    double alpha_prop = prop_alpha[0];

    // Check for numerical stability
    if (!std::isfinite(alpha_prop) || alpha_prop <= 0 || alpha_prop > priorParams[0]) {
      alpha_samples.push_back(alpha_current);
      lambda_samples.push_back(lambda_current);
      continue;
    }

    // Calculate sum(x^alpha) for lambda update
    double sum_x_alpha = 0.0;
    bool valid_sum = true;
    for (int i = 0; i < x.n_rows; i++) {
      if (x(i, 0) > 0) {
        double x_alpha = std::pow(x(i, 0), alpha_prop);
        if (std::isfinite(x_alpha)) {
          sum_x_alpha += x_alpha;
        } else {
          valid_sum = false;
          break;
        }
      }
    }

    if (!valid_sum || sum_x_alpha <= 0) {
      alpha_samples.push_back(alpha_current);
      lambda_samples.push_back(lambda_current);
      continue;
    }

    // Update lambda analytically
    double lambda_prop = 1.0 / R::rgamma(priorParams[1] + x.n_rows,
                                         1.0 / (sum_x_alpha + priorParams[2]));

    // Check lambda validity
    if (!std::isfinite(lambda_prop) || lambda_prop <= 0 || lambda_prop > 1e10) {
      alpha_samples.push_back(alpha_current);
      lambda_samples.push_back(lambda_current);
      continue;
    }

    // Update proposed params with new lambda
    Rcpp::NumericVector prop_lambda(1);
    prop_lambda.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
    prop_lambda[0] = lambda_prop;
    proposed_params[1] = prop_lambda;

    // Calculate proposed likelihood
    double proposed_log_lik = 0.0;
    Rcpp::NumericVector prop_lik_vals = likelihood(arma::vectorise(x), proposed_params);
    for (int k = 0; k < prop_lik_vals.size(); k++) {
      if (prop_lik_vals[k] > 1e-300) {
        proposed_log_lik += std::log(prop_lik_vals[k]);
      } else {
        proposed_log_lik = -std::numeric_limits<double>::infinity();
        break;
      }
    }

    double proposed_log_prior = std::log(priorDensity(proposed_params)[0]);

    // Accept/reject
    double log_ratio = (proposed_log_lik + proposed_log_prior) -
      (current_log_lik + current_log_prior);

    double accept_prob = std::min(1.0, std::exp(log_ratio));
    if (!std::isfinite(accept_prob)) {
      accept_prob = 0.0;
    }

    if (R::runif(0, 1) < accept_prob) {
      current_params = proposed_params;
      alpha_current = alpha_prop;
      lambda_current = lambda_prop;
      current_log_lik = proposed_log_lik;
      current_log_prior = proposed_log_prior;
      accept_count++;
    }

    alpha_samples.push_back(alpha_current);
    lambda_samples.push_back(lambda_current);
  }

  if (iter_count >= MAX_ITER) {
    Rcpp::warning("Weibull posterior draw reached maximum iterations");
  }

  // Return the last n samples
  int actual_samples = alpha_samples.size();
  int start_idx = std::max(0, actual_samples - n);

  Rcpp::NumericVector alpha_final(n);
  Rcpp::NumericVector lambda_final(n);

  for (int i = 0; i < n; i++) {
    int idx = start_idx + i;
    if (idx < actual_samples) {
      alpha_final[i] = alpha_samples[idx];
      lambda_final[i] = lambda_samples[idx];
    } else {
      // Use last available sample
      alpha_final[i] = alpha_samples.back();
      lambda_final[i] = lambda_samples.back();
    }
  }

  alpha_final.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);
  lambda_final.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha_final,
    Rcpp::Named("lambda") = lambda_final
  );
}

void WeibullMixingDistribution::updatePriorParameters(const Rcpp::List& clusterParameters, int n) {
  Rcpp::NumericVector hyperPrior = Rcpp::as<Rcpp::NumericVector>(hyperPriorParameters);
  Rcpp::NumericVector currentPriorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);
  Rcpp::NumericVector alpha_params = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
  Rcpp::NumericVector lambda_params = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

  int numClusters = alpha_params.size();

  // Find max alpha
  double max_alpha = 0.0;
  for (int i = 0; i < numClusters; i++) {
    if (alpha_params[i] > max_alpha) {
      max_alpha = alpha_params[i];
    }
  }

  // Update phi using Pareto distribution
  double xm = std::max(max_alpha, hyperPrior[0]);
  double newPhi = qpareto(R::runif(0, 1), xm, hyperPrior[1] + numClusters);

  // Update gamma parameters
  double sum_inv_lambda = 0.0;
  for (int i = 0; i < numClusters; i++) {
    if (lambda_params[i] > 0) {
      sum_inv_lambda += 1.0 / lambda_params[i];
    }
  }

  double newGamma = R::rgamma(hyperPrior[2] + 2 * numClusters,
                              1.0 / (hyperPrior[3] + sum_inv_lambda));

  Rcpp::NumericMatrix newPriorParams(1, 3);
  newPriorParams(0, 0) = newPhi;
  newPriorParams(0, 1) = currentPriorParams[1];
  newPriorParams(0, 2) = newGamma;

  priorParameters = newPriorParams;
}

// Helper function for Pareto quantile
double WeibullMixingDistribution::qpareto(double p, double xm, double alpha) const {
  if (p < 0 || p > 1) return R_NaN;
  return xm * std::pow(1 - p, -1.0 / alpha);
}

// NonConjugateWeibullDP implementation
NonConjugateWeibullDP::NonConjugateWeibullDP() : mixingDistribution(nullptr), numberClusters(0), m(3) {
  // Constructor
}

NonConjugateWeibullDP::~NonConjugateWeibullDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void NonConjugateWeibullDP::clusterComponentUpdate() {
  int n = data.n_rows;

  for (int i = 0; i < n; i++) {
    int currentLabel = clusterLabels[i];

    // Validate current label
    if (currentLabel < 0 || currentLabel >= (int)pointsPerCluster.n_elem) {
      Rcpp::Rcerr << "Warning: Invalid cluster label " << currentLabel
                  << " for point " << i << ". Resetting to 0.\n";
      currentLabel = 0;
      clusterLabels[i] = 0;
    }

    // Create auxiliary parameters
    Rcpp::List aux;
    bool needAux = (pointsPerCluster[currentLabel] == 1);

    if (needAux) {
      // Generate m auxiliary parameters
      aux = mixingDistribution->priorDraw(m);
    } else {
      // Generate m-1 auxiliary parameters
      aux = mixingDistribution->priorDraw(m - 1);
    }

    // Calculate probabilities for each possible cluster (including auxiliary)
    int totalOptions = numberClusters + m - 1;
    if (!needAux) totalOptions = numberClusters + m - 1;

    Rcpp::NumericVector probs(totalOptions);

    // Existing clusters
    for (int j = 0; j < numberClusters; j++) {
      double count = (j == currentLabel) ?
      pointsPerCluster[j] - 1 : pointsPerCluster[j];

      if (count > 0 || j != currentLabel) {
        Rcpp::NumericVector alpha_j = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
        Rcpp::NumericVector lambda_j = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

        Rcpp::NumericVector alpha_temp(1), lambda_temp(1);
        alpha_temp.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        lambda_temp.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

        if (j < alpha_j.size() && j < lambda_j.size()) {
          alpha_temp[0] = alpha_j[j];
          lambda_temp[0] = lambda_j[j];

          Rcpp::List theta_j = Rcpp::List::create(alpha_temp, lambda_temp);
          arma::vec xi = data.row(i).t();
          Rcpp::NumericVector lik = mixingDistribution->likelihood(xi, theta_j);

          probs[j] = count * lik[0];
        } else {
          probs[j] = 0.0;
        }
      } else {
        probs[j] = 0.0;
      }
    }

    // Auxiliary clusters
    Rcpp::NumericVector aux_alpha = aux[0];
    Rcpp::NumericVector aux_lambda = aux[1];

    for (int j = 0; j < m - 1; j++) {
      int prob_idx = numberClusters + j;
      if (prob_idx < probs.size() && j < aux_alpha.size()) {
        Rcpp::NumericVector alpha_temp(1), lambda_temp(1);
        alpha_temp.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        lambda_temp.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        alpha_temp[0] = aux_alpha[j];
        lambda_temp[0] = aux_lambda[j];

        Rcpp::List theta_aux = Rcpp::List::create(alpha_temp, lambda_temp);
        arma::vec xi = data.row(i).t();
        Rcpp::NumericVector lik = mixingDistribution->likelihood(xi, theta_aux);

        probs[prob_idx] = (alpha / (m - 1)) * lik[0];
      }
    }

    // Handle numerical issues
    double prob_sum = 0.0;
    for (int j = 0; j < probs.size(); j++) {
      if (!std::isfinite(probs[j]) || probs[j] < 0) {
        probs[j] = 0.0;
      }
      prob_sum += probs[j];
    }

    if (prob_sum <= 0) {
      // Fallback to uniform
      for (int j = 0; j < probs.size(); j++) {
        probs[j] = 1.0 / probs.size();
      }
      prob_sum = 1.0;
    }

    // Normalize
    probs = probs / prob_sum;

    // Sample new label
    int newLabel = 0;
    double u = R::runif(0, 1);
    double cumProb = 0.0;
    for (int j = 0; j < probs.size(); j++) {
      cumProb += probs[j];
      if (u <= cumProb) {
        newLabel = j;
        break;
      }
    }

    // Update the state
    pointsPerCluster[currentLabel]--;

    // Perform the label change
    Rcpp::List updateResult = clusterLabelChange(i, newLabel, currentLabel, aux);

    // Update state from result
    clusterLabels = Rcpp::as<arma::uvec>(updateResult["clusterLabels"]);
    pointsPerCluster = Rcpp::as<arma::uvec>(updateResult["pointsPerCluster"]);
    clusterParameters = updateResult["clusterParameters"];
    numberClusters = updateResult["numberClusters"];
  }
}

void NonConjugateWeibullDP::clusterParameterUpdate() {
  for (int k = 0; k < numberClusters; k++) {
    arma::uvec clusterIndices = arma::find(clusterLabels == k);
    if (clusterIndices.n_elem > 0) {
      arma::mat clusterData = data.rows(clusterIndices);

      // Extract current parameters
      Rcpp::NumericVector alpha_params = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector lambda_params = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      if (k < alpha_params.size() && k < lambda_params.size()) {
        // Draw from posterior
        Rcpp::List postDraw = mixingDistribution->posteriorDraw(clusterData, 1);

        Rcpp::NumericVector new_alpha = postDraw[0];
        Rcpp::NumericVector new_lambda = postDraw[1];

        alpha_params[k] = new_alpha[0];
        lambda_params[k] = new_lambda[0];

        clusterParameters[0] = alpha_params;
        clusterParameters[1] = lambda_params;
      }
    }
  }
}

void NonConjugateWeibullDP::updateAlpha() {
  // Same implementation as Beta/Normal
  double x = R::rbeta(alpha + 1.0, n);
  Rcpp::NumericVector currentAlphaPrior = Rcpp::as<Rcpp::NumericVector>(alphaPriorParameters);

  double log_x = std::log(x);
  double pi1 = currentAlphaPrior[0] + numberClusters - 1.0;
  double pi2 = n * (currentAlphaPrior[1] - log_x);

  double pi_val = pi1 / (pi1 + pi2);
  if (!std::isfinite(pi_val)) {
    pi_val = 0.5;
  }

  double postShape;
  if (R::runif(0, 1) < pi_val) {
    postShape = currentAlphaPrior[0] + numberClusters;
  } else {
    postShape = currentAlphaPrior[0] + numberClusters - 1.0;
  }

  double postRate = currentAlphaPrior[1] - log_x;
  if (postRate <= 0) postRate = 1e-6;

  alpha = R::rgamma(postShape, 1.0 / postRate);
  if (alpha <= 0) alpha = 1e-6;
}

Rcpp::List NonConjugateWeibullDP::clusterLabelChange(int i, int newLabel, int currentLabel,
                                                     const Rcpp::List& aux) {
  if (newLabel == currentLabel) {
    // No change needed
    pointsPerCluster[currentLabel]++;
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  // Extract current parameters
  Rcpp::NumericVector alpha_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
  Rcpp::NumericVector lambda_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

  // Handle new cluster assignment
  if (newLabel < numberClusters) {
    // Existing cluster
    pointsPerCluster[newLabel]++;
    clusterLabels[i] = newLabel;
  } else {
    // New cluster from auxiliary
    int aux_idx = newLabel - numberClusters;
    Rcpp::NumericVector aux_alpha = aux[0];
    Rcpp::NumericVector aux_lambda = aux[1];

    if (aux_idx < aux_alpha.size()) {
      // Add new cluster
      numberClusters++;

      // Extend parameter vectors
      Rcpp::NumericVector new_alpha_vec(numberClusters);
      Rcpp::NumericVector new_lambda_vec(numberClusters);

      for (int j = 0; j < alpha_vec.size(); j++) {
        new_alpha_vec[j] = alpha_vec[j];
        new_lambda_vec[j] = lambda_vec[j];
      }

      new_alpha_vec[numberClusters - 1] = aux_alpha[aux_idx];
      new_lambda_vec[numberClusters - 1] = aux_lambda[aux_idx];

      alpha_vec = new_alpha_vec;
      lambda_vec = new_lambda_vec;

      // Extend pointsPerCluster
      arma::uvec new_points = arma::zeros<arma::uvec>(numberClusters);
      for (int j = 0; j < (int)pointsPerCluster.n_elem; j++) {
        new_points[j] = pointsPerCluster[j];
      }
      new_points[numberClusters - 1] = 1;
      pointsPerCluster = new_points;

      clusterLabels[i] = numberClusters - 1;
    }
  }

  // Clean up empty clusters
  if (pointsPerCluster[currentLabel] == 0 && numberClusters > 1) {
    // Remove empty cluster
    numberClusters--;

    // Create new vectors without the empty cluster
    Rcpp::NumericVector new_alpha_vec(numberClusters);
    Rcpp::NumericVector new_lambda_vec(numberClusters);
    arma::uvec new_points = arma::zeros<arma::uvec>(numberClusters);

    int new_idx = 0;
    std::map<int, int> label_map;

    for (int j = 0; j < alpha_vec.size(); j++) {
      if (j != currentLabel && j < (int)pointsPerCluster.n_elem) {
        new_alpha_vec[new_idx] = alpha_vec[j];
        new_lambda_vec[new_idx] = lambda_vec[j];
        new_points[new_idx] = pointsPerCluster[j];
        label_map[j] = new_idx;
        new_idx++;
      }
    }

    // Update all cluster labels
    for (int k = 0; k < (int)clusterLabels.size(); k++) {
      if (label_map.find(clusterLabels[k]) != label_map.end()) {
        clusterLabels[k] = label_map[clusterLabels[k]];
      }
    }

    alpha_vec = new_alpha_vec;
    lambda_vec = new_lambda_vec;
    pointsPerCluster = new_points;
  }

  // Update cluster parameters
  clusterParameters = Rcpp::List::create(alpha_vec, lambda_vec);

  return Rcpp::List::create(
    Rcpp::Named("clusterLabels") = clusterLabels,
    Rcpp::Named("pointsPerCluster") = pointsPerCluster,
    Rcpp::Named("clusterParameters") = clusterParameters,
    Rcpp::Named("numberClusters") = numberClusters
  );
}

} // namespace dp
