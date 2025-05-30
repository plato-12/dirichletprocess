// src/MVNormal2Distribution.cpp
#include "../inst/include/MVNormal2Distribution.h"
#include "../inst/include/RcppConversions.h"
#include <RcppArmadillo.h>
#include <cmath>

namespace dp {

// MVNormal2MixingDistribution implementation
MVNormal2MixingDistribution::MVNormal2MixingDistribution(const Rcpp::List& priorParams) {
  distribution = "mvnormal2";
  conjugate = false;
  priorParameters = priorParams;

  // Extract prior parameters
  if (priorParams.containsElementNamed("mu0")) {
    SEXP mu0_sexp = priorParams["mu0"];
    if (Rf_isMatrix(mu0_sexp)) {
      arma::mat temp_mu0_mat = Rcpp::as<arma::mat>(mu0_sexp); // Convert to arma::mat
      if (temp_mu0_mat.n_rows == 1) { // If R matrix is 1xN (already a row vector shape)
        mu0 = temp_mu0_mat; // Assign directly (arma::mat to arma::rowvec if mat is 1xN)
      } else if (temp_mu0_mat.n_cols == 1) { // If R matrix is Nx1 (a column vector shape)
        mu0 = temp_mu0_mat.t(); // Transpose to 1xN and assign
      } else {
        Rcpp::stop("mu0 in priorParams, if a matrix, must be a row or column vector.");
      }
    } else { // It's an R vector (NumericVector)
      // Rcpp::as<arma::vec> converts an R vector to an Armadillo column vector
      arma::vec temp_mu0_col_vec = Rcpp::as<arma::vec>(mu0_sexp);
      mu0 = temp_mu0_col_vec.t(); // Transpose the column vector to a row vector for mu0
    }
  }

  if (priorParams.containsElementNamed("sigma0")) {
    sigma0 = Rcpp::as<arma::mat>(priorParams["sigma0"]);
  }

  if (priorParams.containsElementNamed("phi0")) {
    phi0 = Rcpp::as<arma::mat>(priorParams["phi0"]);
  }

  if (priorParams.containsElementNamed("nu0")) {
    nu0 = Rcpp::as<double>(priorParams["nu0"]);
  }

  // Set default MH step size if not provided
  if (!priorParameters.containsElementNamed("mhStepSize")) {
    mhStepSize = Rcpp::NumericVector::create(1.0, 1.0);
  }
}

MVNormal2MixingDistribution::~MVNormal2MixingDistribution() {
  // Destructor
}

Rcpp::NumericVector MVNormal2MixingDistribution::likelihood(const arma::vec& x, const Rcpp::List& theta) const {
  // Extract parameters from theta
  Rcpp::NumericVector mu_array = theta[0];
  Rcpp::NumericVector sig_array = theta[1];

  // Get dimensions
  Rcpp::IntegerVector mu_dim = mu_array.attr("dim");
  Rcpp::IntegerVector sig_dim = sig_array.attr("dim");

  int d = mu_dim[1];  // Number of dimensions
  int n_clusters = mu_dim[2];  // Number of clusters

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

    // Calculate multivariate normal likelihood
    double log_det_val;
    double sign;
    arma::log_det(log_det_val, sign, sig_k);

    if (sign <= 0) {
      result[k] = 1e-300;
      continue;
    }

    arma::mat sig_inv;
    try {
      sig_inv = arma::inv_sympd(sig_k);
    } catch(...) {
      result[k] = 1e-300;
      continue;
    }

    double log_const = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det_val;
    arma::vec x_centered = x - mu_k;
    double quad_form = arma::as_scalar(x_centered.t() * sig_inv * x_centered);
    result[k] = std::exp(log_const - 0.5 * quad_form);
  }

  return result;
}

Rcpp::List MVNormal2MixingDistribution::priorDraw(int n) const {
  int d = mu0.n_elem;

  Rcpp::NumericVector mu_arr = Rcpp::NumericVector(Rcpp::Dimension(1, d, n));
  Rcpp::NumericVector sig_arr = Rcpp::NumericVector(Rcpp::Dimension(d, d, n));

  for (int i = 0; i < n; i++) {
    // Draw Sigma from Inverse-Wishart (corrected parameterization)
    arma::mat sig_draw = arma::iwishrnd(phi0, nu0);

    // Draw mu from Multivariate Normal given Sigma
    arma::vec mu_draw = arma::mvnrnd(mu0.t(), sigma0);

    // Store in arrays
    for (int j = 0; j < d; j++) {
      mu_arr[j + i * d] = mu_draw(j);
    }

    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        sig_arr[j + k * d + i * d * d] = sig_draw(j, k);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sig") = sig_arr
  );
}

Rcpp::List MVNormal2MixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  if (!x.is_finite()) {
    Rcpp::stop("Input data contains non-finite values");
  }

  int d = x.n_cols;
  if (d == 0 || x.n_rows == 0) {
    // Return prior draw if no data
    return priorDraw(n);
  }

  // Arrays to store results
  Rcpp::NumericVector mu_arr = Rcpp::NumericVector(Rcpp::Dimension(1, d, n));
  Rcpp::NumericVector sig_arr = Rcpp::NumericVector(Rcpp::Dimension(d, d, n));

  // Initialize with a reasonable starting value
  arma::vec mu_samp = arma::mean(x, 0).t();

  for (int i = 0; i < n; i++) {
    // Update Sigma given current mu
    double nu_n = x.n_rows + nu0;
    arma::mat phi_n = phi0;

    for (arma::uword j = 0; j < x.n_rows; j++) {
      arma::vec diff = x.row(j).t() - mu_samp;
      phi_n += diff * diff.t();
    }

    // Draw new Sigma
    arma::mat sig_samp = arma::iwishrnd(arma::inv_sympd(phi_n), nu_n);

    // Update mu given new Sigma
    arma::mat sig_n = arma::inv_sympd(arma::inv_sympd(sigma0) + x.n_rows * arma::inv_sympd(sig_samp));
    arma::vec mu_n = sig_n * (x.n_rows * arma::inv_sympd(sig_samp) * arma::mean(x, 0).t() +
      arma::inv_sympd(sigma0) * mu0.t());

    // Draw new mu
    mu_samp = arma::mvnrnd(mu_n, sig_n);

    // Store results
    for (int j = 0; j < d; j++) {
      mu_arr[j + i * d] = mu_samp(j);
    }

    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        sig_arr[j + k * d + i * d * d] = sig_samp(j, k);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sig") = sig_arr
  );
}

// NonConjugateMVNormal2DP implementation
NonConjugateMVNormal2DP::NonConjugateMVNormal2DP() : mixingDistribution(nullptr), numberClusters(0), m(3) {
  // Constructor
}

NonConjugateMVNormal2DP::~NonConjugateMVNormal2DP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void NonConjugateMVNormal2DP::clusterComponentUpdate() {
  int n = data.n_rows;

  for (int i = 0; i < n; i++) {
    int currentLabel = clusterLabels[i];

    // Remove point from current cluster
    pointsPerCluster[currentLabel]--;

    // Generate auxiliary parameters
    Rcpp::List aux;
    if (pointsPerCluster[currentLabel] == 0) {
      // If cluster is now empty, we need m-1 auxiliary parameters
      aux = mixingDistribution->priorDraw(m - 1);

      // Include the current cluster's parameters as one of the auxiliary
      Rcpp::NumericVector mu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector sig_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      Rcpp::NumericVector mu_aux = aux[0];
      Rcpp::NumericVector sig_aux = aux[1];

      // Get dimensions
      Rcpp::IntegerVector mu_dim = mu_vec.attr("dim");
      int d = mu_dim[1];

      // Create new arrays including current cluster params
      Rcpp::NumericVector mu_combined = Rcpp::NumericVector(Rcpp::Dimension(1, d, m));
      Rcpp::NumericVector sig_combined = Rcpp::NumericVector(Rcpp::Dimension(d, d, m));

      // Copy current cluster parameters
      for (int j = 0; j < d; j++) {
        mu_combined[j] = mu_vec[j + currentLabel * d];
      }
      for (int j = 0; j < d; j++) {
        for (int k = 0; k < d; k++) {
          sig_combined[j + k * d] = sig_vec[j + k * d + currentLabel * d * d];
        }
      }

      // Copy auxiliary parameters
      for (int idx = 1; idx < m; idx++) {
        for (int j = 0; j < d; j++) {
          mu_combined[j + idx * d] = mu_aux[j + (idx-1) * d];
        }
        for (int j = 0; j < d; j++) {
          for (int k = 0; k < d; k++) {
            sig_combined[j + k * d + idx * d * d] = sig_aux[j + k * d + (idx-1) * d * d];
          }
        }
      }

      aux = Rcpp::List::create(mu_combined, sig_combined);
    } else {
      // Generate m auxiliary parameters
      aux = mixingDistribution->priorDraw(m);
    }

    // Calculate probabilities
    int totalLabels = numberClusters + m;
    Rcpp::NumericVector probs(totalLabels);

    // Existing clusters
    for (int j = 0; j < numberClusters; j++) {
      if (pointsPerCluster[j] > 0) {
        // Create parameter list for cluster j
        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::Named("mu") = clusterParameters[0],
                                               Rcpp::Named("sig") = clusterParameters[1]
        );

        Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), clusterParam);
        probs[j] = pointsPerCluster[j] * lik[j];
      } else {
        probs[j] = 0.0;
      }
    }

    // Auxiliary clusters
    for (int j = 0; j < m; j++) {
      Rcpp::List auxParam = Rcpp::List::create(
        Rcpp::Named("mu") = aux[0],
                               Rcpp::Named("sig") = aux[1]
      );

      Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), auxParam);
      probs[numberClusters + j] = (alpha / m) * lik[j];
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

    // Update cluster assignment
    Rcpp::List updateResult = clusterLabelChange(i, newLabel, currentLabel, aux);

    // Update state from result
    clusterLabels = Rcpp::as<arma::uvec>(updateResult["clusterLabels"]);
    pointsPerCluster = Rcpp::as<arma::uvec>(updateResult["pointsPerCluster"]);
    clusterParameters = updateResult["clusterParameters"];
    numberClusters = updateResult["numberClusters"];
  }
}

void NonConjugateMVNormal2DP::clusterParameterUpdate() {
  for (int k = 0; k < numberClusters; k++) {
    arma::uvec clusterIndices = arma::find(clusterLabels == k);
    if (clusterIndices.n_elem > 0) {
      arma::mat clusterData = data.rows(clusterIndices);

      // Draw from posterior
      Rcpp::List postDraw = mixingDistribution->posteriorDraw(clusterData, mhDraws);

      // Update cluster parameters - extract last sample
      Rcpp::NumericVector mu_samples = postDraw[0];
      Rcpp::NumericVector sig_samples = postDraw[1];

      // Get dimensions
      Rcpp::IntegerVector mu_dim = mu_samples.attr("dim");
      int d = mu_dim[1];

      // Extract current parameter arrays
      Rcpp::NumericVector mu_params = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector sig_params = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      // Update with last sample
      int last_idx = mhDraws - 1;
      for (int j = 0; j < d; j++) {
        mu_params[j + k * d] = mu_samples[j + last_idx * d];
      }
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          sig_params[i + j * d + k * d * d] = sig_samples[i + j * d + last_idx * d * d];
        }
      }

      clusterParameters[0] = mu_params;
      clusterParameters[1] = sig_params;
    }
  }
}

void NonConjugateMVNormal2DP::updateAlpha() {
  // Same implementation as other non-conjugate cases
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

Rcpp::List NonConjugateMVNormal2DP::clusterLabelChange(int i, int newLabel, int currentLabel,
                                                       const Rcpp::List& aux) {
  if (newLabel == currentLabel) {
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  // Extract current parameters
  Rcpp::NumericVector mu_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
  Rcpp::NumericVector sig_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

  // Get dimensions
  Rcpp::IntegerVector mu_dim = mu_vec.attr("dim");
  int d = mu_dim[1];

  // 1. Remove point from old cluster (already done by caller)

  // 2. Assign to new cluster
  if (newLabel < numberClusters) {
    // Existing cluster
    pointsPerCluster[newLabel]++;
    clusterLabels[i] = newLabel;

    // If old cluster is now empty, remove it
    if (pointsPerCluster[currentLabel] == 0) {
      numberClusters--;
      pointsPerCluster.shed_row(currentLabel);

      // Create new arrays with reduced size
      Rcpp::NumericVector new_mu = Rcpp::NumericVector(Rcpp::Dimension(1, d, numberClusters));
      Rcpp::NumericVector new_sig = Rcpp::NumericVector(Rcpp::Dimension(d, d, numberClusters));

      // Copy parameters, skipping the removed cluster
      int new_k = 0;
      for (int k = 0; k < numberClusters + 1; k++) {
        if (k != currentLabel) {
          for (int j = 0; j < d; j++) {
            new_mu[j + new_k * d] = mu_vec[j + k * d];
          }
          for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
              new_sig[i + j * d + new_k * d * d] = sig_vec[i + j * d + k * d * d];
            }
          }
          new_k++;
        }
      }

      clusterParameters[0] = new_mu;
      clusterParameters[1] = new_sig;

      // Update labels
      for (arma::uword j = 0; j < clusterLabels.n_elem; j++) {
        if (clusterLabels[j] > (unsigned int)currentLabel) {
          clusterLabels[j]--;
        }
      }
    }
  } else {
    // New cluster from auxiliary parameters
    int auxIndex = newLabel - numberClusters;

    if (pointsPerCluster[currentLabel] == 0) {
      // Replace empty cluster with auxiliary
      Rcpp::NumericVector aux_mu = aux[0];
      Rcpp::NumericVector aux_sig = aux[1];

      // Copy auxiliary parameters to current cluster position
      for (int j = 0; j < d; j++) {
        mu_vec[j + currentLabel * d] = aux_mu[j + auxIndex * d];
      }
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          sig_vec[i + j * d + currentLabel * d * d] =
            aux_sig[i + j * d + auxIndex * d * d];
        }
      }

      pointsPerCluster[currentLabel] = 1;
      clusterLabels[i] = currentLabel;
    } else {
      // Create new cluster
      Rcpp::NumericVector aux_mu = aux[0];
      Rcpp::NumericVector aux_sig = aux[1];

      // Create expanded arrays
      Rcpp::NumericVector new_mu = Rcpp::NumericVector(Rcpp::Dimension(1, d, numberClusters + 1));
      Rcpp::NumericVector new_sig = Rcpp::NumericVector(Rcpp::Dimension(d, d, numberClusters + 1));

      // Copy existing parameters
      for (int k = 0; k < numberClusters; k++) {
        for (int j = 0; j < d; j++) {
          new_mu[j + k * d] = mu_vec[j + k * d];
        }
        for (int i = 0; i < d; i++) {
          for (int j = 0; j < d; j++) {
            new_sig[i + j * d + k * d * d] = sig_vec[i + j * d + k * d * d];
          }
        }
      }

      // Add new cluster parameters
      for (int j = 0; j < d; j++) {
        new_mu[j + numberClusters * d] = aux_mu[j + auxIndex * d];
      }
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          new_sig[i + j * d + numberClusters * d * d] =
            aux_sig[i + j * d + auxIndex * d * d];
        }
      }

      clusterParameters[0] = new_mu;
      clusterParameters[1] = new_sig;

      clusterLabels[i] = numberClusters;
      pointsPerCluster.resize(numberClusters + 1);
      pointsPerCluster[numberClusters] = 1;
      numberClusters++;
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
