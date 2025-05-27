// src/ExponentialDistribution.cpp
#include "../inst/include/ExponentialDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <cmath>

namespace dp {

// ExponentialMixingDistribution implementation
ExponentialMixingDistribution::ExponentialMixingDistribution(const Rcpp::NumericVector& priorParams) {
  distribution = "exponential";
  conjugate = true;
  priorParameters = priorParams;
}

ExponentialMixingDistribution::~ExponentialMixingDistribution() {
  // Destructor
}

Rcpp::NumericVector ExponentialMixingDistribution::likelihood(const arma::vec& x, const Rcpp::List& theta) const {
  // Extract rate parameter from theta
  Rcpp::NumericVector lambda_array = theta[0];

  int n_data = x.n_elem;
  Rcpp::NumericVector result(n_data);

  // Get the first element (assuming single cluster for now)
  double lambda = lambda_array[0];

  if (lambda <= 0) {
    result.fill(1e-300);
    return result;
  }

  // Calculate exponential likelihood: f(x|λ) = λ * exp(-λ * x)
  for (int i = 0; i < n_data; i++) {
    if (x[i] >= 0) {
      result[i] = lambda * std::exp(-lambda * x[i]);
    } else {
      result[i] = 0.0; // Exponential distribution is only defined for x >= 0
    }
  }

  return result;
}

Rcpp::List ExponentialMixingDistribution::priorDraw(int n) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  // Prior parameters: alpha0 (shape), beta0 (rate)
  double alpha0 = priorParams[0];
  double beta0 = priorParams[1];

  Rcpp::NumericVector lambda(n);

  // Draw from Gamma(alpha0, beta0)
  for (int i = 0; i < n; i++) {
    lambda[i] = R::rgamma(alpha0, 1.0/beta0); // Note: R::rgamma uses scale = 1/rate
  }

  // Convert to 3D array with dimension (1,1,n)
  Rcpp::NumericVector lambda_arr(n);
  lambda_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    lambda_arr[i] = lambda[i];
  }

  return Rcpp::List::create(Rcpp::Named("lambda") = lambda_arr);
}

Rcpp::NumericMatrix ExponentialMixingDistribution::posteriorParameters(const arma::mat& x) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  int n_x = x.n_rows;
  double sum_x = arma::sum(arma::vectorise(x));

  double alpha0 = priorParams[0];
  double beta0 = priorParams[1];

  // Posterior parameters for Gamma distribution
  double alpha_n = alpha0 + n_x;
  double beta_n = beta0 + sum_x;

  Rcpp::NumericMatrix result(1, 2);
  result(0, 0) = alpha_n;
  result(0, 1) = beta_n;

  return result;
}

Rcpp::List ExponentialMixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  // First compute posterior parameters
  Rcpp::NumericMatrix postParams = posteriorParameters(x);

  double alpha_n = postParams(0, 0);
  double beta_n = postParams(0, 1);

  Rcpp::NumericVector lambda(n);

  // Draw from posterior Gamma(alpha_n, beta_n)
  for (int i = 0; i < n; i++) {
    lambda[i] = R::rgamma(alpha_n, 1.0/beta_n);
  }

  // Convert to 3D array with dimension (1,1,n)
  Rcpp::NumericVector lambda_arr(n);
  lambda_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    lambda_arr[i] = lambda[i];
  }

  return Rcpp::List::create(Rcpp::Named("lambda") = lambda_arr);
}

Rcpp::NumericVector ExponentialMixingDistribution::predictive(const arma::vec& x) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  int n = x.n_elem;
  Rcpp::NumericVector result(n);

  double alpha0 = priorParams[0];
  double beta0 = priorParams[1];

  for (int i = 0; i < n; i++) {
    // For single observation
    double alpha_post = alpha0 + 1;
    double beta_post = beta0 + x[i];

    // Predictive distribution calculation
    result[i] = (R::gammafn(alpha_post) / R::gammafn(alpha0)) *
      std::pow(beta0, alpha0) / std::pow(beta_post, alpha_post);
  }

  return result;
}

// ConjugateExponentialDP implementation
ConjugateExponentialDP::ConjugateExponentialDP(Rcpp::List dpObj) {
  // Initialize from the R list object
  data = Rcpp::as<arma::mat>(dpObj["data"]);
  n = data.n_rows;
  alpha = dpObj["alpha"];
  clusterLabels = Rcpp::as<arma::uvec>(dpObj["clusterLabels"]);
  pointsPerCluster = Rcpp::as<arma::uvec>(dpObj["pointsPerCluster"]);
  numberClusters = dpObj["numberClusters"];
  clusterParameters = dpObj["clusterParameters"];
  predictiveArray = Rcpp::as<arma::vec>(dpObj["predictiveArray"]);

  Rcpp::List mixingDistributionList = dpObj["mixingDistribution"];
  Rcpp::NumericVector priorParams = mixingDistributionList["priorParameters"];
  mixingDistribution = new ExponentialMixingDistribution(priorParams);
}


ConjugateExponentialDP::~ConjugateExponentialDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void ConjugateExponentialDP::initialisePredictive() {
  // Calculate predictive probabilities for all data points
  predictiveArray = mixingDistribution->predictive(arma::vectorise(data));
}

void ConjugateExponentialDP::clusterComponentUpdate() {
  // Implementation of Chinese Restaurant Process for conjugate case
  int n = data.n_rows;

  for (int i = 0; i < n; i++) {
    int currentLabel = clusterLabels[i];

    // Remove point from current cluster
    pointsPerCluster[currentLabel]--;

    // Calculate probabilities for existing clusters
    Rcpp::NumericVector probs(numberClusters + 1);

    // Probability for existing clusters
    for (int j = 0; j < numberClusters; j++) {
      if (pointsPerCluster[j] > 0) {
        // Extract parameters for cluster j
        Rcpp::NumericVector lambda_vec = clusterParameters[0];

        // Create properly formatted parameter array
        Rcpp::NumericVector lambda_j(1);
        lambda_j[0] = lambda_vec[j];
        lambda_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);

        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::Named("lambda") = lambda_j
        );

        double likelihood = mixingDistribution->likelihood(data.row(i).t(), clusterParam)[0];
        probs[j] = pointsPerCluster[j] * likelihood;
      } else {
        probs[j] = 0.0;
      }
    }

    // Probability for new cluster
    probs[numberClusters] = alpha * predictiveArray[i];

    // Normalize probabilities
    double probSum = arma::sum(arma::vec(probs));
    if (probSum == 0) {
      // If all probabilities are 0, make them uniform
      probs.fill(1.0 / probs.size());
    } else {
      probs = probs / probSum;
    }

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

    // Restore point count before the change
    pointsPerCluster[currentLabel]++;

    // Update cluster assignment using clusterLabelChange
    Rcpp::List updateResult = clusterLabelChange(i, newLabel, currentLabel);

    // Update state from result
    clusterLabels = Rcpp::as<arma::uvec>(updateResult["clusterLabels"]);
    pointsPerCluster = Rcpp::as<arma::uvec>(updateResult["pointsPerCluster"]);
    clusterParameters = updateResult["clusterParameters"];
    numberClusters = updateResult["numberClusters"];
  }
}

void ConjugateExponentialDP::clusterParameterUpdate() {
  // Update parameters for each cluster
  for (int k = 0; k < numberClusters; k++) {
    // Get data points assigned to this cluster
    arma::uvec clusterIndices = arma::find(clusterLabels == k);

    if (clusterIndices.n_elem > 0) {
      arma::mat clusterData = data.rows(clusterIndices);

      // Draw from posterior
      Rcpp::List postDraw = mixingDistribution->posteriorDraw(clusterData, 1);

      // Update cluster parameters
      Rcpp::NumericVector lambda_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector new_lambda = postDraw["lambda"];

      lambda_vec[k] = new_lambda[0];
      clusterParameters[0] = lambda_vec;
    }
  }
}

void ConjugateExponentialDP::updateAlpha() {
  // Implementation of alpha update using auxiliary variable method
  double x = R::rbeta(alpha + 1.0, n);

  // Cast alphaPriorParameters to NumericVector
  Rcpp::NumericVector alphaPriors = Rcpp::as<Rcpp::NumericVector>(alphaPriorParameters);

  double pi1 = alphaPriors[0] + numberClusters - 1.0;
  double pi2 = n * (alphaPriors[1] - log(x));
  double pi_ratio = pi1 / (pi1 + pi2);

  double postShape, postRate;
  if (R::runif(0, 1) < pi_ratio) {
    postShape = alphaPriors[0] + numberClusters;
  } else {
    postShape = alphaPriors[0] + numberClusters - 1.0;
  }
  postRate = alphaPriors[1] - log(x);

  alpha = R::rgamma(postShape, 1.0/postRate);
}

Rcpp::List ConjugateExponentialDP::clusterLabelChange(int i, int newLabel, int currentLabel) {
  if (newLabel == currentLabel) {
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  arma::mat x_i = data.row(i);

  // 1. Remove point from old cluster
  pointsPerCluster[currentLabel]--;

  // 2. Assign point to new cluster
  clusterLabels[i] = newLabel;
  if (newLabel == numberClusters) { // This is a new cluster
    numberClusters++;
    pointsPerCluster.resize(numberClusters);
    pointsPerCluster(newLabel) = 1;

    // Safely create copy of parameter vector, modify, and assign back
    Rcpp::NumericVector lambda_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));

    Rcpp::List postDraw = mixingDistribution->posteriorDraw(x_i, 1);
    lambda_vec.push_back(Rcpp::as<Rcpp::NumericVector>(postDraw["lambda"])[0]);

    clusterParameters[0] = lambda_vec;

  } else { // This is an existing cluster
    pointsPerCluster[newLabel]++;
  }

  // 3. If the old cluster is now empty, remove it
  if (pointsPerCluster[currentLabel] == 0) {
    pointsPerCluster.shed_row(currentLabel);

    // Safely create copy, modify, and assign back
    Rcpp::NumericVector lambda_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
    lambda_vec.erase(currentLabel);
    clusterParameters[0] = lambda_vec;

    numberClusters--;

    // Shift all labels that were greater than the removed cluster's label
    for (arma::uword j = 0; j < clusterLabels.n_elem; j++) {
      if (clusterLabels[j] > (unsigned int)currentLabel) {
        clusterLabels[j]--;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("clusterLabels") = clusterLabels,
    Rcpp::Named("pointsPerCluster") = pointsPerCluster,
    Rcpp::Named("clusterParameters") = clusterParameters,
    Rcpp::Named("numberClusters") = numberClusters
  );
}

} // namespace dp
