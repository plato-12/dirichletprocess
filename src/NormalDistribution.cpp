// src/NormalDistribution.cpp
#include "../inst/include/NormalDistribution.h"
#include "../inst/include/RcppConversions.h"

namespace dp {

// NormalMixingDistribution implementation
NormalMixingDistribution::NormalMixingDistribution(const Rcpp::NumericVector& priorParams) {
  distribution = "normal";
  conjugate = true;
  priorParameters = priorParams;
}

NormalMixingDistribution::~NormalMixingDistribution() {
  // Destructor
}

Rcpp::NumericVector NormalMixingDistribution::likelihood(const arma::vec& x, const Rcpp::List& theta) const {
  // Extract parameters from theta
  Rcpp::NumericVector mu_array = theta[0];
  Rcpp::NumericVector sigma_array = theta[1];

  int n_clusters = mu_array.size();
  int n_data = x.n_elem;

  Rcpp::NumericVector result(n_data);

  // For now, use first cluster only
  if (n_clusters > 0 && n_data > 0) {
    double mu = mu_array[0];
    double sigma = sigma_array[0];

    for (int i = 0; i < n_data; i++) {
      result[i] = R::dnorm(x[i], mu, sigma, false);
    }
  }

  return result;
}

Rcpp::List NormalMixingDistribution::priorDraw(int n) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  // Prior parameters: mu0, kappa0, alpha0, beta0
  double mu0 = priorParams[0];
  double kappa0 = priorParams[1];
  double alpha0 = priorParams[2];
  double beta0 = priorParams[3];

  Rcpp::NumericVector mu(n);
  Rcpp::NumericVector sigma(n);
  Rcpp::NumericVector lambda_vec(n);

  // 1. Draw all n lambda values first to match R's vectorized behavior
  for (int i = 0; i < n; i++) {
    lambda_vec[i] = R::rgamma(alpha0, 1.0/beta0);
  }

  // 2. Then draw all n mu values
  for (int i = 0; i < n; i++) {
    double lambda = lambda_vec[i];
    // Draw mu from Normal(mu0, 1/(kappa0*lambda))
    double mu_sd = 1.0/sqrt(kappa0 * lambda);
    mu[i] = R::rnorm(mu0, mu_sd);

    // sigma = sqrt(1/lambda)
    sigma[i] = sqrt(1.0/lambda);
  }

  // Convert to 3D arrays with dimension (1,1,n)
  Rcpp::NumericVector mu_arr(n);
  Rcpp::NumericVector sigma_arr(n);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);
  sigma_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    mu_arr[i] = mu[i];
    sigma_arr[i] = sigma[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sigma") = sigma_arr
  );
}

Rcpp::List NormalMixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  // First compute posterior parameters
  Rcpp::NumericMatrix postParams = posteriorParameters(x);

  Rcpp::NumericVector mu(n);
  Rcpp::NumericVector sigma(n);
  Rcpp::NumericVector lambda_vec(n);

  // Extract posterior parameters
  double mu_n = postParams(0, 0);
  double kappa_n = postParams(0, 1);
  double alpha_n = postParams(0, 2);
  double beta_n = postParams(0, 3);

  // 1. Draw all n lambda values first
  for (int i = 0; i < n; i++) {
    lambda_vec[i] = R::rgamma(alpha_n, 1.0/beta_n);
  }

  // 2. Then draw all n mu values
  for (int i = 0; i < n; i++) {
    double lambda = lambda_vec[i];
    // Draw mu from Normal(mu_n, 1/(kappa_n*lambda))
    double mu_sd = 1.0/sqrt(kappa_n * lambda);
    mu[i] = R::rnorm(mu_n, mu_sd);

    // sigma = sqrt(1/lambda)
    sigma[i] = sqrt(1.0/lambda);
  }

  // Convert to 3D arrays with dimension (1,1,n)
  Rcpp::NumericVector mu_arr(n);
  Rcpp::NumericVector sigma_arr(n);
  mu_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);
  sigma_arr.attr("dim") = Rcpp::IntegerVector::create(1, 1, n);

  for (int i = 0; i < n; i++) {
    mu_arr[i] = mu[i];
    sigma_arr[i] = sigma[i];
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sigma") = sigma_arr
  );
}

Rcpp::NumericMatrix NormalMixingDistribution::posteriorParameters(const arma::mat& x) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  int n_x = x.n_rows;
  double ybar = arma::mean(arma::vectorise(x));

  double mu0 = priorParams[0];
  double kappa0 = priorParams[1];
  double alpha0 = priorParams[2];
  double beta0 = priorParams[3];

  double mu_n = (kappa0 * mu0 + n_x * ybar) / (kappa0 + n_x);
  double kappa_n = kappa0 + n_x;
  double alpha_n = alpha0 + n_x / 2.0;
  double beta_n = beta0 + 0.5 * arma::sum(arma::square(arma::vectorise(x) - ybar)) +
    kappa0 * n_x * std::pow(ybar - mu0, 2) / (2.0 * (kappa0 + n_x));

  Rcpp::NumericMatrix result(1, 4);
  result(0, 0) = mu_n;
  result(0, 1) = kappa_n;
  result(0, 2) = alpha_n;
  result(0, 3) = beta_n;

  return result;
}

Rcpp::NumericVector NormalMixingDistribution::predictive(const arma::vec& x) const {
  Rcpp::NumericVector priorParams = Rcpp::as<Rcpp::NumericVector>(priorParameters);

  int n = x.n_elem;
  Rcpp::NumericVector result(n);

  for (int i = 0; i < n; i++) {
    Rcpp::NumericMatrix postParams = posteriorParameters(arma::mat(&x[i], 1, 1));

    double predictive_val = (R::gammafn(postParams(0, 2)) / R::gammafn(priorParams[2])) *
      (std::pow(priorParams[3], priorParams[2]) /
        std::pow(postParams(0, 3), postParams(0, 2))) *
          std::sqrt(priorParams[1] / postParams(0, 1));

    result[i] = predictive_val;
  }

  return result;
}

// Static methods for direct testing
Rcpp::List NormalMixingDistribution::priorDrawStatic(const Rcpp::NumericVector& priorParams, int n) {
  NormalMixingDistribution md(priorParams);
  return md.priorDraw(n);
}

Rcpp::List NormalMixingDistribution::posteriorDrawStatic(const Rcpp::NumericVector& priorParams, const arma::mat& x, int n) {
  NormalMixingDistribution md(priorParams);
  return md.posteriorDraw(x, n);
}

// ConjugateNormalDP implementation
ConjugateNormalDP::ConjugateNormalDP() : mixingDistribution(nullptr), numberClusters(0) {
  // Constructor
}

ConjugateNormalDP::~ConjugateNormalDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void ConjugateNormalDP::clusterComponentUpdate() {
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
        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::as<Rcpp::NumericVector>(clusterParameters[0])[j],
                                                             Rcpp::as<Rcpp::NumericVector>(clusterParameters[1])[j]
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

    // Put the point's count back before the change.
    // The actual state change is handled entirely by clusterLabelChange now.
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

void ConjugateNormalDP::clusterParameterUpdate() {
  // Update parameters for each cluster
  for (int k = 0; k < numberClusters; k++) {
    // Get data points assigned to this cluster
    arma::uvec clusterIndices = arma::find(clusterLabels == k);

    if (clusterIndices.n_elem > 0) {
      arma::mat clusterData = data.rows(clusterIndices);

      // Draw from posterior
      Rcpp::List postDraw = mixingDistribution->posteriorDraw(clusterData, 1);

      // Update cluster parameters
      Rcpp::NumericVector mu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector sigma_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      Rcpp::NumericVector new_mu = postDraw["mu"];
      Rcpp::NumericVector new_sigma = postDraw["sigma"];

      mu_vec[k] = new_mu[0];
      sigma_vec[k] = new_sigma[0];

      clusterParameters[0] = mu_vec;
      clusterParameters[1] = sigma_vec;
    }
  }
}

void ConjugateNormalDP::updateAlpha() {
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

Rcpp::List ConjugateNormalDP::clusterLabelChange(int i, int newLabel, int currentLabel) {
  if (newLabel == currentLabel) {
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  arma::mat x_i = data.row(i);

  // 1. Remove point from its old cluster
  pointsPerCluster[currentLabel]--;

  // 2. Assign point to its new cluster and update parameters
  clusterLabels[i] = newLabel;
  if (newLabel == numberClusters) { // This is a new cluster
    numberClusters++;
    pointsPerCluster.resize(numberClusters);
    pointsPerCluster(newLabel) = 1;

    // Safely create copies of parameter vectors, modify, and assign back
    Rcpp::NumericVector mu_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
    Rcpp::NumericVector sigma_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

    Rcpp::List postDraw = mixingDistribution->posteriorDraw(x_i, 1);
    mu_vec.push_back(Rcpp::as<Rcpp::NumericVector>(postDraw["mu"])[0]);
    sigma_vec.push_back(Rcpp::as<Rcpp::NumericVector>(postDraw["sigma"])[0]);

    clusterParameters[0] = mu_vec;
    clusterParameters[1] = sigma_vec;

  } else { // This is an existing cluster
    pointsPerCluster[newLabel]++;
  }

  // 3. If the old cluster is now empty, remove it
  if (pointsPerCluster[currentLabel] == 0) {
    pointsPerCluster.shed_row(currentLabel);

    // Safely create copies, modify, and assign back
    Rcpp::NumericVector mu_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
    Rcpp::NumericVector sigma_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

    mu_vec.erase(currentLabel);
    sigma_vec.erase(currentLabel);

    clusterParameters[0] = mu_vec;
    clusterParameters[1] = sigma_vec;

    numberClusters--;

    // Shift all labels that were greater than the removed cluster's label
    for (arma::uword j = 0; j < clusterLabels.n_elem; j++) {
      if (clusterLabels[j] > currentLabel) {
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


void ConjugateNormalDP::initialisePredictive() {
  // Calculate predictive probabilities for all data points
  predictiveArray = mixingDistribution->predictive(arma::vectorise(data));
}

} // namespace dp
