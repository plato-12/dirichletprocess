// src/WeibullDistribution.cpp
#include "../inst/include/WeibullDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <cmath>
#include <limits>

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
    if (x[i] >= 0 && lambda > 0 && alpha > 0 && !std::isinf(lambda)) {
      double y = std::pow(lambda, -1.0) * alpha * std::pow(x[i], alpha - 1.0) *
        std::exp(-std::pow(lambda, -1.0) * std::pow(x[i], alpha));
      result[i] = (y > 0 && std::isfinite(y)) ? y : 1e-300;
    } else {
      result[i] = (x[i] < 0) ? 0.0 : 1e-300;
    }
  }

  return result;
}

Rcpp::List WeibullMixingDistribution::priorDraw(int n) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  Rcpp::NumericVector alpha_values(n);
  Rcpp::NumericVector lambda_values(n);

  // Fix: priorParams indexing (0-based in C++)
  // priorParams[0] = phi, priorParams[1] = alpha0, priorParams[2] = beta0
  for (int i = 0; i < n; i++) {
    alpha_values[i] = R::runif(0.0, priorParams[0]);
    // R code: lambdas <- 1/rgamma(n, priorParameters[2], priorParameters[3])
    // R's rgamma uses shape and rate, but R::rgamma uses shape and scale
    // Need to convert rate to scale: scale = 1/rate
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
  double new_alpha = std::abs(alpha_old + mhStep[0] * R::rnorm(0.0, 1.7)); // Match R implementation

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
  // Use the special Weibull Metropolis-Hastings with analytical lambda update
  Rcpp::List start_pos = priorDraw(1);

  int mhDraws = 100; // Default number of MH iterations
  if (priorParameters.hasAttribute("mhDraws")) {
    mhDraws = Rcpp::as<int>(priorParameters.attr("mhDraws"));
  }

  // Perform Metropolis-Hastings sampling
  Rcpp::NumericVector alpha_samples(mhDraws);
  Rcpp::NumericVector lambda_samples(mhDraws);

  Rcpp::NumericVector current_alpha = start_pos[0];
  double alpha_current = current_alpha[0];

  // Initialize with analytical lambda update
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);
  double sum_x_alpha = 0.0;
  for (arma::uword j = 0; j < x.n_rows; j++) {
    sum_x_alpha += std::pow(x(j, 0), alpha_current);
  }
  double lambda_current = 1.0 / R::rgamma(x.n_rows + priorParams[1],
                                          1.0 / (sum_x_alpha + priorParams[2]));

  alpha_samples[0] = alpha_current;
  lambda_samples[0] = lambda_current;

  // Current likelihood and prior
  Rcpp::NumericVector current_params_alpha(1);
  Rcpp::NumericVector current_params_lambda(1);
  current_params_alpha.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  current_params_lambda.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
  current_params_alpha[0] = alpha_current;
  current_params_lambda[0] = lambda_current;
  Rcpp::List current_params = Rcpp::List::create(current_params_alpha, current_params_lambda);

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

  // Metropolis-Hastings loop
  for (int iter = 1; iter < mhDraws; iter++) {
    // Propose new alpha
    Rcpp::List proposed_params = mhParameterProposal(current_params);
    Rcpp::NumericVector prop_alpha = proposed_params[0];
    double alpha_prop = prop_alpha[0];

    // Analytically update lambda given proposed alpha
    sum_x_alpha = 0.0;
    for (arma::uword j = 0; j < x.n_rows; j++) {
      sum_x_alpha += std::pow(x(j, 0), alpha_prop);
    }
    double lambda_prop = 1.0 / R::rgamma(x.n_rows + priorParams[1],
                                         1.0 / (sum_x_alpha + priorParams[2]));

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

    alpha_samples[iter] = alpha_current;
    lambda_samples[iter] = lambda_current;
  }

  // Return the last n samples
  Rcpp::NumericVector alpha_final(n);
  Rcpp::NumericVector lambda_final(n);

  for (int i = 0; i < n; i++) {
    int idx = mhDraws - n + i;
    if (idx >= 0) {
      alpha_final[i] = alpha_samples[idx];
      lambda_final[i] = lambda_samples[idx];
    } else {
      alpha_final[i] = alpha_samples[mhDraws - 1];
      lambda_final[i] = lambda_samples[mhDraws - 1];
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

    // Ensure currentLabel is a valid index for pointsPerCluster
    if (currentLabel < 0 || currentLabel >= (int)pointsPerCluster.n_elem) {
      Rcpp::stop("Invalid cluster label encountered for point %d: %d (max allowed: %d)", i, currentLabel, pointsPerCluster.n_elem - 1);
    }

    // Create a copy of the current state for probability calculations
    arma::uvec tempPointsPerCluster = pointsPerCluster;
    tempPointsPerCluster[currentLabel]--;

    // Generate auxiliary parameters
    Rcpp::List aux;
    if (tempPointsPerCluster[currentLabel] == 0) {
      // If cluster would be empty, we need m-1 auxiliary parameters
      aux = mixingDistribution->priorDraw(m - 1);

      // Include the current cluster's parameters as one of the auxiliary
      Rcpp::NumericVector alpha_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector lambda_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      Rcpp::NumericVector alpha_aux = aux[0];
      Rcpp::NumericVector lambda_aux = aux[1];

      // Create new arrays including current cluster params
      Rcpp::NumericVector alpha_combined(m);
      Rcpp::NumericVector lambda_combined(m);
      alpha_combined.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);
      lambda_combined.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);

      alpha_combined[0] = alpha_vec[currentLabel];
      lambda_combined[0] = lambda_vec[currentLabel];

      for (int j = 0; j < m-1; j++) {
        alpha_combined[j+1] = alpha_aux[j];
        lambda_combined[j+1] = lambda_aux[j];
      }

      aux = Rcpp::List::create(alpha_combined, lambda_combined);
    } else {
      // Generate m auxiliary parameters
      aux = mixingDistribution->priorDraw(m);
    }

    // Calculate probabilities using temporary counts
    int totalLabels = numberClusters + m;
    Rcpp::NumericVector probs(totalLabels);

    // Existing clusters
    for (int j = 0; j < numberClusters; j++) {
      if (tempPointsPerCluster[j] > 0) {
        // Extract parameters for cluster j
        Rcpp::NumericVector alpha_vec = clusterParameters[0];
        Rcpp::NumericVector lambda_vec = clusterParameters[1];

        // Create properly formatted parameter arrays
        Rcpp::NumericVector alpha_j(1);
        Rcpp::NumericVector lambda_j(1);
        alpha_j[0] = alpha_vec[j];
        lambda_j[0] = lambda_vec[j];

        // Add dimension attributes
        alpha_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        lambda_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::Named("alpha") = alpha_j,
          Rcpp::Named("lambda") = lambda_j
        );

        Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), clusterParam);
        probs[j] = tempPointsPerCluster[j] * lik[0]; // Use temp counts
      } else {
        probs[j] = 0.0;
      }
    }

    // Auxiliary clusters
    for (int j = 0; j < m; j++) {
      Rcpp::NumericVector alpha_aux = aux[0];
      Rcpp::NumericVector lambda_aux = aux[1];

      Rcpp::NumericVector alpha_j(1);
      Rcpp::NumericVector lambda_j(1);
      alpha_j[0] = alpha_aux[j];
      lambda_j[0] = lambda_aux[j];

      alpha_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      lambda_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

      Rcpp::List auxParam = Rcpp::List::create(
        Rcpp::Named("alpha") = alpha_j,
        Rcpp::Named("lambda") = lambda_j
      );

      Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), auxParam);
      probs[numberClusters + j] = (alpha / m) * lik[0];
    }

    // Handle edge cases
    if (Rcpp::is_true(Rcpp::any(Rcpp::is_nan(probs)))) {
      for (int j = 0; j < probs.size(); j++) {
        if (std::isnan(probs[j])) probs[j] = 0.0;
      }
    }

    if (Rcpp::is_true(Rcpp::all(probs == 0))) {
      probs.fill(1.0 / probs.size());
    }

    // Normalize
    double probSum = Rcpp::sum(probs);
    probs = probs / probSum;

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

    // Now update the actual state
    // First decrement the count from current cluster
    pointsPerCluster[currentLabel]--;

    // Then perform the label change
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

      Rcpp::NumericVector alpha_start(1);
      Rcpp::NumericVector lambda_start(1);
      alpha_start[0] = alpha_params[k];
      lambda_start[0] = lambda_params[k];
      alpha_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      lambda_start.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

      Rcpp::List start_pos = Rcpp::List::create(
        Rcpp::Named("alpha") = alpha_start,
        Rcpp::Named("lambda") = lambda_start
      );

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
    // No change needed, but still need to re-increment the count
    pointsPerCluster[currentLabel]++;
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  // Note: pointsPerCluster[currentLabel] has already been decremented

  // Extract current parameters
  Rcpp::NumericVector alpha_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
  Rcpp::NumericVector lambda_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

  // Assign to new cluster
  if (newLabel < numberClusters) {
    // Existing cluster
    pointsPerCluster[newLabel]++;
    clusterLabels[i] = newLabel;

    // If old cluster is now empty, remove it
    if (pointsPerCluster[currentLabel] == 0 && currentLabel != newLabel) {
      numberClusters--;

      // Remove the empty cluster from pointsPerCluster
      arma::uvec new_points = arma::uvec(numberClusters);
      int idx = 0;
      for (int j = 0; j < (int)pointsPerCluster.n_elem; j++) {
        if (j != currentLabel) {
          new_points[idx++] = pointsPerCluster[j];
        }
      }
      pointsPerCluster = new_points;

      // Remove from parameter vectors
      Rcpp::NumericVector new_alpha_vec(numberClusters);
      Rcpp::NumericVector new_lambda_vec(numberClusters);

      idx = 0;
      for (int j = 0; j < alpha_vec.size(); j++) {
        if (j != currentLabel) {
          new_alpha_vec[idx] = alpha_vec[j];
          new_lambda_vec[idx] = lambda_vec[j];
          idx++;
        }
      }

      alpha_vec = new_alpha_vec;
      lambda_vec = new_lambda_vec;

      // Update all labels that were greater than currentLabel
      for (arma::uword j = 0; j < clusterLabels.n_elem; j++) {
        if ((int)clusterLabels[j] > currentLabel) {
          clusterLabels[j]--;
        }
      }

      // If the point was assigned to a label that got shifted, update it
      if (newLabel > currentLabel) {
        clusterLabels[i] = newLabel - 1;
      }
    }
  } else {
    // New cluster from auxiliary parameters
    int auxIndex = newLabel - numberClusters;

    if (pointsPerCluster[currentLabel] == 0) {
      // Replace empty cluster with auxiliary
      Rcpp::NumericVector aux_alpha = aux[0];
      Rcpp::NumericVector aux_lambda = aux[1];

      alpha_vec[currentLabel] = aux_alpha[auxIndex];
      lambda_vec[currentLabel] = aux_lambda[auxIndex];

      pointsPerCluster[currentLabel] = 1;
      clusterLabels[i] = currentLabel;
    } else {
      // Create new cluster
      Rcpp::NumericVector aux_alpha = aux[0];
      Rcpp::NumericVector aux_lambda = aux[1];

      // Expand vectors
      Rcpp::NumericVector new_alpha_vec(numberClusters + 1);
      Rcpp::NumericVector new_lambda_vec(numberClusters + 1);

      for (int j = 0; j < numberClusters; j++) {
        new_alpha_vec[j] = alpha_vec[j];
        new_lambda_vec[j] = lambda_vec[j];
      }

      new_alpha_vec[numberClusters] = aux_alpha[auxIndex];
      new_lambda_vec[numberClusters] = aux_lambda[auxIndex];

      alpha_vec = new_alpha_vec;
      lambda_vec = new_lambda_vec;

      // Expand pointsPerCluster
      arma::uvec new_points = arma::uvec(numberClusters + 1);
      for (int j = 0; j < numberClusters; j++) {
        new_points[j] = pointsPerCluster[j];
      }
      new_points[numberClusters] = 1;
      pointsPerCluster = new_points;

      clusterLabels[i] = numberClusters;
      numberClusters++;
    }
  }

  // Update cluster parameters
  clusterParameters[0] = alpha_vec;
  clusterParameters[1] = lambda_vec;

  return Rcpp::List::create(
    Rcpp::Named("clusterLabels") = clusterLabels,
    Rcpp::Named("pointsPerCluster") = pointsPerCluster,
    Rcpp::Named("clusterParameters") = clusterParameters,
    Rcpp::Named("numberClusters") = numberClusters
  );
}

Rcpp::List NonConjugateWeibullDP::metropolisHastings(const arma::mat& x, const Rcpp::List& startPos, int noDraws) {
  return mixingDistribution->posteriorDraw(x, noDraws);
}

} // namespace dp
