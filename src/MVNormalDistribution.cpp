// src/MVNormalDistribution.cpp
#include "../inst/include/MVNormalDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <RcppArmadillo.h>

namespace dp {

// MVNormalMixingDistribution implementation
MVNormalMixingDistribution::MVNormalMixingDistribution(const Rcpp::List& priorParams) {
  distribution = "mvnormal";
  conjugate = true;
  priorParameters = priorParams;

  // Extract prior parameters
  if (priorParams.containsElementNamed("mu0")) {
    Rcpp::NumericVector mu0_vec = Rcpp::as<Rcpp::NumericVector>(priorParams["mu0"]);
    mu0 = arma::vec(mu0_vec.begin(), mu0_vec.size());
  }

  if (priorParams.containsElementNamed("kappa0")) {
    kappa0 = Rcpp::as<double>(priorParams["kappa0"]);
  }

  if (priorParams.containsElementNamed("Lambda")) {
    Lambda = Rcpp::as<arma::mat>(priorParams["Lambda"]);
  }

  if (priorParams.containsElementNamed("nu")) {
    nu = Rcpp::as<double>(priorParams["nu"]);
  }
}

MVNormalMixingDistribution::~MVNormalMixingDistribution() {
  // Destructor
}

arma::vec MVNormalMixingDistribution::mvnLikelihood(const arma::mat& x,
                                                    const arma::vec& mu,
                                                    const arma::mat& sigma) const {
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec result(n);

  // Calculate log-determinant and inverse once
  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, sigma);

  if (sign <= 0) {
    // Sigma is not positive definite
    result.fill(1e-300);
    return result;
  }

  arma::mat sigma_inv;
  try {
    sigma_inv = arma::inv_sympd(sigma);
  } catch(...) {
    // If inversion fails, return very small likelihood
    result.fill(1e-300);
    return result;
  }

  double log_const = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det_val;

  for (int i = 0; i < n; i++) {
    arma::vec x_centered = x.row(i).t() - mu;
    double quad_form = arma::as_scalar(x_centered.t() * sigma_inv * x_centered);
    result(i) = std::exp(log_const - 0.5 * quad_form);
  }

  return result;
}

Rcpp::NumericVector MVNormalMixingDistribution::likelihood(const arma::vec& x,
                                                           const Rcpp::List& theta) const {
  // Extract parameters - handle the array structure
  Rcpp::NumericVector mu_array = theta["mu"];
  Rcpp::NumericVector sig_array = theta["sig"];

  // Get dimensions
  Rcpp::IntegerVector mu_dim = mu_array.attr("dim");
  Rcpp::IntegerVector sig_dim = sig_array.attr("dim");

  int d = mu_dim[1]; // Number of dimensions
  int n_clusters = mu_dim[2]; // Number of clusters

  // Convert x to matrix (single row)
  arma::mat x_mat(1, x.n_elem);
  x_mat.row(0) = x.t();

  Rcpp::NumericVector result(n_clusters);

  for (int k = 0; k < n_clusters; k++) {
    // Extract mu for cluster k
    arma::vec mu_k(d);
    for (int j = 0; j < d; j++) {
      mu_k(j) = mu_array[j + k * d];
    }

    // Extract sigma for cluster k
    arma::mat sig_k(d, d);
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        sig_k(i, j) = sig_array[i + j * d + k * d * d];
      }
    }

    arma::vec lik = mvnLikelihood(x_mat, mu_k, sig_k);
    result[k] = lik(0);
  }

  return result;
}

Rcpp::List MVNormalMixingDistribution::priorDraw(int n) const {
  int d = mu0.n_elem;

  // Arrays to store results
  Rcpp::NumericVector mu_arr = Rcpp::NumericVector(Rcpp::Dimension(1, d, n));
  Rcpp::NumericVector sig_arr = Rcpp::NumericVector(Rcpp::Dimension(d, d, n));

  for (int i = 0; i < n; i++) {
    // Draw precision matrix from Wishart (to match R implementation)
    arma::mat prec_draw = arma::wishrnd(Lambda, nu);

    // Convert to covariance matrix (Sigma = Lambda^{-1})
    arma::mat sig_draw = arma::inv_sympd(prec_draw);

    // Draw mu from Multivariate Normal given Sigma
    arma::mat cov_mu = sig_draw / kappa0;
    arma::vec mu_draw = arma::mvnrnd(mu0, cov_mu);

    // Store in arrays
    for (int j = 0; j < d; j++) {
      mu_arr[j + i * d] = mu_draw(j);
    }

    // Store the precision matrix (not the covariance) to match R implementation
    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        sig_arr[j + k * d + i * d * d] = prec_draw(j, k);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sig") = sig_arr
  );
}

Rcpp::List MVNormalMixingDistribution::posteriorParameters(const arma::mat& x) const {
  int n = x.n_rows;
  int d = x.n_cols;

  // Handle empty data case
  if (n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("mu_n") = mu0,
      Rcpp::Named("t_n") = Lambda,  // Changed from Lambda_n to t_n
      Rcpp::Named("kappa_n") = kappa0,
      Rcpp::Named("nu_n") = nu
    );
  }

  // Compute sample statistics
  arma::vec x_bar = arma::mean(x, 0).t();

  // Posterior parameters
  double kappa_n = kappa0 + n;
  arma::vec mu_n = (kappa0 * mu0 + n * x_bar) / kappa_n;
  double nu_n = nu + n;

  // Compute scatter matrix
  arma::mat S = arma::zeros(d, d);
  if (n > 1) {
    S = (n - 1) * arma::cov(x);
  }

  // Update Lambda (called t_n in R code)
  arma::vec diff = x_bar - mu0;
  arma::mat t_n = Lambda + S + (kappa0 * n / kappa_n) * (diff * diff.t());

  return Rcpp::List::create(
    Rcpp::Named("mu_n") = mu_n,
    Rcpp::Named("t_n") = t_n,  // Changed from Lambda_n to t_n
    Rcpp::Named("kappa_n") = kappa_n,
    Rcpp::Named("nu_n") = nu_n
  );
}

Rcpp::List MVNormalMixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  // Get posterior parameters
  Rcpp::List post_params = posteriorParameters(x);

  arma::vec mu_n = Rcpp::as<arma::vec>(post_params["mu_n"]);
  arma::mat t_n = Rcpp::as<arma::mat>(post_params["t_n"]);  // Changed from Lambda_n
  double kappa_n = Rcpp::as<double>(post_params["kappa_n"]);
  double nu_n = Rcpp::as<double>(post_params["nu_n"]);

  int d = mu_n.n_elem;

  // Arrays to store results
  Rcpp::NumericVector mu_arr = Rcpp::NumericVector(Rcpp::Dimension(1, d, n));
  Rcpp::NumericVector sig_arr = Rcpp::NumericVector(Rcpp::Dimension(d, d, n));

  for (int i = 0; i < n; i++) {
    // Draw precision from Wishart (to match R implementation)
    arma::mat prec_draw = arma::wishrnd(t_n, nu_n);

    // Convert to covariance
    arma::mat sig_draw = arma::inv_sympd(prec_draw);

    // Draw mu from Multivariate Normal given Sigma
    arma::mat cov_mu = sig_draw / kappa_n;
    arma::vec mu_draw = arma::mvnrnd(mu_n, cov_mu);

    // Store in arrays
    for (int j = 0; j < d; j++) {
      mu_arr[j + i * d] = mu_draw(j);
    }

    // Store precision matrix to match R
    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        sig_arr[j + k * d + i * d * d] = prec_draw(j, k);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sig") = sig_arr
  );
}

Rcpp::NumericVector MVNormalMixingDistribution::predictive(const arma::mat& x) const {
  int n = x.n_rows;
  int d = x.n_cols;
  Rcpp::NumericVector result(n);

  double pi_const = std::pow(M_PI, -0.5 * d);

  for (int i = 0; i < n; i++) {
    arma::mat x_i = x.row(i);
    Rcpp::List post_params = posteriorParameters(x_i);

    arma::vec mu_n = Rcpp::as<arma::vec>(post_params["mu_n"]);
    arma::mat t_n = Rcpp::as<arma::mat>(post_params["t_n"]);  // Changed from Lambda_n
    double kappa_n = Rcpp::as<double>(post_params["kappa_n"]);
    double nu_n = Rcpp::as<double>(post_params["nu_n"]);

    double det_Lambda, det_t_n;
    double sign;
    arma::log_det(det_Lambda, sign, Lambda);
    det_Lambda = std::exp(det_Lambda) * sign;
    arma::log_det(det_t_n, sign, t_n);
    det_t_n = std::exp(det_t_n) * sign;

    double ratio_det = std::pow(det_Lambda / det_t_n, nu / 2.0);
    double ratio_kappa = std::pow(kappa0 / kappa_n, d / 2.0);

    // Compute multivariate gamma ratio
    double gamma_ratio = 1.0;
    for (int j = 0; j < d; j++) {
      gamma_ratio *= R::gammafn((nu_n + 1.0 - j) / 2.0) / R::gammafn((nu + 1.0 - j) / 2.0);
    }

    result[i] = pi_const * ratio_kappa * ratio_det * gamma_ratio;
  }

  return result;
}

// Static methods
Rcpp::List MVNormalMixingDistribution::priorDrawStatic(const Rcpp::List& priorParams, int n) {
  MVNormalMixingDistribution md(priorParams);
  return md.priorDraw(n);
}

Rcpp::List MVNormalMixingDistribution::posteriorDrawStatic(const Rcpp::List& priorParams,
                                                           const arma::mat& x, int n) {
  MVNormalMixingDistribution md(priorParams);
  return md.posteriorDraw(x, n);
}

// ConjugateMVNormalDP implementation
ConjugateMVNormalDP::ConjugateMVNormalDP() : mixingDistribution(nullptr), numberClusters(0) {
  // Constructor
}

ConjugateMVNormalDP::~ConjugateMVNormalDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void ConjugateMVNormalDP::initialisePredictive() {
  // Calculate predictive probabilities for all data points
  predictiveArray = mixingDistribution->predictive(data);
}

void ConjugateMVNormalDP::clusterComponentUpdate() {
  // Implementation similar to Normal but for multivariate case
  int n = data.n_rows;

  for (int i = 0; i < n; i++) {
    int currentLabel = clusterLabels[i];

    // Remove point from current cluster
    pointsPerCluster[currentLabel]--;

    // Calculate probabilities for existing clusters
    Rcpp::NumericVector probs(numberClusters + 1);

    // Get parameters from clusterParameters list
    Rcpp::NumericVector mu_array = clusterParameters["mu"];
    Rcpp::NumericVector sig_array = clusterParameters["sig"];

    // Probability for existing clusters
    for (int j = 0; j < numberClusters; j++) {
      if (pointsPerCluster[j] > 0) {
        // Create parameter list for cluster j
        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::Named("mu") = mu_array,
          Rcpp::Named("sig") = sig_array
        );

        Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), clusterParam);
        probs[j] = pointsPerCluster[j] * lik[j];
      } else {
        probs[j] = 0.0;
      }
    }

    // Probability for new cluster
    probs[numberClusters] = alpha * predictiveArray[i];

    // Handle edge cases
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

    // Restore point count before calling clusterLabelChange
    pointsPerCluster[currentLabel]++;

    // Update cluster assignment
    Rcpp::List updateResult = clusterLabelChange(i, newLabel, currentLabel);

    // Update state from result
    clusterLabels = Rcpp::as<arma::uvec>(updateResult["clusterLabels"]);
    pointsPerCluster = Rcpp::as<arma::uvec>(updateResult["pointsPerCluster"]);
    clusterParameters = updateResult["clusterParameters"];
    numberClusters = updateResult["numberClusters"];
  }
}

void ConjugateMVNormalDP::clusterParameterUpdate() {
  // Update parameters for each cluster
  for (int k = 0; k < numberClusters; k++) {
    // Get data points assigned to this cluster
    arma::uvec clusterIndices = arma::find(clusterLabels == k);

    if (clusterIndices.n_elem > 0) {
      arma::mat clusterData = data.rows(clusterIndices);

      // Draw from posterior
      Rcpp::List postDraw = mixingDistribution->posteriorDraw(clusterData, 1);

      // Update cluster parameters - this is more complex for multivariate case
      // Need to handle the array structure properly
      Rcpp::NumericVector mu_array = clusterParameters["mu"];
      Rcpp::NumericVector sig_array = clusterParameters["sig"];

      Rcpp::NumericVector new_mu = postDraw["mu"];
      Rcpp::NumericVector new_sig = postDraw["sig"];

      // Get dimensions
      Rcpp::IntegerVector mu_dim = mu_array.attr("dim");
      int d = mu_dim[1];

      // Update the k-th cluster parameters
      for (int j = 0; j < d; j++) {
        mu_array[j + k * d] = new_mu[j];
      }

      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          sig_array[i + j * d + k * d * d] = new_sig[i + j * d];
        }
      }

      clusterParameters["mu"] = mu_array;
      clusterParameters["sig"] = sig_array;
    }
  }
}

void ConjugateMVNormalDP::updateAlpha() {
  // Same implementation as Normal case
  double x = R::rbeta(alpha + 1.0, n);

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

Rcpp::List ConjugateMVNormalDP::clusterLabelChange(int i, int newLabel, int currentLabel) {
  if (newLabel == currentLabel) {
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  arma::mat x_i = data.row(i);

  // Extract current parameters
  Rcpp::NumericVector mu_array = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters["mu"]));
  Rcpp::NumericVector sig_array = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters["sig"]));

  // Get dimensions
  Rcpp::IntegerVector mu_dim = mu_array.attr("dim");
  int d = mu_dim[1];

  // 1. Remove point from old cluster
  pointsPerCluster[currentLabel]--;

  // 2. Assign point to new cluster
  clusterLabels[i] = newLabel;

  if (newLabel == numberClusters) {
    // New cluster
    numberClusters++;
    pointsPerCluster.resize(numberClusters);
    pointsPerCluster(newLabel) = 1;

    // Draw parameters for new cluster
    Rcpp::List postDraw = mixingDistribution->posteriorDraw(x_i, 1);
    Rcpp::NumericVector new_mu = postDraw["mu"];
    Rcpp::NumericVector new_sig = postDraw["sig"];

    // Create new arrays with expanded size
    Rcpp::NumericVector new_mu_array = Rcpp::NumericVector(Rcpp::Dimension(1, d, numberClusters));
    Rcpp::NumericVector new_sig_array = Rcpp::NumericVector(Rcpp::Dimension(d, d, numberClusters));

    // Copy existing parameters
    for (int k = 0; k < numberClusters - 1; k++) {
      for (int j = 0; j < d; j++) {
        new_mu_array[j + k * d] = mu_array[j + k * d];
      }
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          new_sig_array[i + j * d + k * d * d] = sig_array[i + j * d + k * d * d];
        }
      }
    }

    // Add new cluster parameters
    for (int j = 0; j < d; j++) {
      new_mu_array[j + (numberClusters - 1) * d] = new_mu[j];
    }
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        new_sig_array[i + j * d + (numberClusters - 1) * d * d] = new_sig[i + j * d];
      }
    }

    clusterParameters["mu"] = new_mu_array;
    clusterParameters["sig"] = new_sig_array;

  } else {
    // Existing cluster
    pointsPerCluster[newLabel]++;
  }

  // 3. If old cluster is empty, remove it
  if (pointsPerCluster[currentLabel] == 0) {
    pointsPerCluster.shed_row(currentLabel);
    numberClusters--;

    // Create new arrays with reduced size
    Rcpp::NumericVector new_mu_array = Rcpp::NumericVector(Rcpp::Dimension(1, d, numberClusters));
    Rcpp::NumericVector new_sig_array = Rcpp::NumericVector(Rcpp::Dimension(d, d, numberClusters));

    // Copy parameters, skipping the removed cluster
    int new_k = 0;
    for (int k = 0; k < numberClusters + 1; k++) {
      if (k != currentLabel) {
        for (int j = 0; j < d; j++) {
          new_mu_array[j + new_k * d] = mu_array[j + k * d];
        }
        for (int i = 0; i < d; i++) {
          for (int j = 0; j < d; j++) {
            new_sig_array[i + j * d + new_k * d * d] = sig_array[i + j * d + k * d * d];
          }
        }
        new_k++;
      }
    }

    clusterParameters["mu"] = new_mu_array;
    clusterParameters["sig"] = new_sig_array;

    // Shift labels
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


