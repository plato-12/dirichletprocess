// src/BetaDP.cpp
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"

namespace dp {

// NonConjugateBetaDP constructor and destructor
NonConjugateBetaDP::NonConjugateBetaDP() : mixingDistribution(nullptr), numberClusters(0), m(3) {
  // Constructor
}

NonConjugateBetaDP::~NonConjugateBetaDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
    mixingDistribution = nullptr;
  }
}

void NonConjugateBetaDP::clusterComponentUpdate() {
  int n = data.n_rows;

  for (int i = 0; i < n; i++) {
    int currentLabel = clusterLabels[i];

    // Validate currentLabel
    if (currentLabel >= numberClusters) {
      Rcpp::stop("Invalid cluster label encountered");
    }

    // Remove point from current cluster temporarily
    pointsPerCluster[currentLabel]--;

    // Generate auxiliary parameters
    Rcpp::List aux;
    bool currentClusterEmpty = (pointsPerCluster[currentLabel] == 0);

    if (currentClusterEmpty) {
      // Current cluster is empty, include it as auxiliary
      aux = mixingDistribution->priorDraw(m - 1);

      // Include current cluster params as first auxiliary
      Rcpp::NumericVector mu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]);
      Rcpp::NumericVector nu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      Rcpp::NumericVector mu_aux = aux[0];
      Rcpp::NumericVector nu_aux = aux[1];

      Rcpp::NumericVector mu_combined(m);
      Rcpp::NumericVector nu_combined(m);

      mu_combined[0] = mu_vec[currentLabel];
      nu_combined[0] = nu_vec[currentLabel];

      for (int j = 1; j < m; j++) {
        mu_combined[j] = mu_aux[j-1];
        nu_combined[j] = nu_aux[j-1];
      }

      mu_combined.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);
      nu_combined.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);

      aux = Rcpp::List::create(mu_combined, nu_combined);
    } else {
      // Generate m auxiliary parameters
      aux = mixingDistribution->priorDraw(m);
    }

    // Calculate probabilities
    Rcpp::NumericVector probs(numberClusters + m);

    // Existing clusters
    for (int j = 0; j < numberClusters; j++) {
      if (pointsPerCluster[j] > 0) {
        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::NumericVector::create(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0])[j]),
          Rcpp::NumericVector::create(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1])[j])
        );

        // Set dimensions
        Rcpp::NumericVector mu_j = clusterParam[0];
        Rcpp::NumericVector nu_j = clusterParam[1];
        mu_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        nu_j.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
        clusterParam[0] = mu_j;
        clusterParam[1] = nu_j;

        Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), clusterParam);
        probs[j] = pointsPerCluster[j] * lik[0];
      } else {
        probs[j] = 0.0;
      }
    }

    // Auxiliary clusters
    for (int j = 0; j < m; j++) {
      Rcpp::List auxParam = Rcpp::List::create(
        Rcpp::NumericVector::create(Rcpp::as<Rcpp::NumericVector>(aux[0])[j]),
        Rcpp::NumericVector::create(Rcpp::as<Rcpp::NumericVector>(aux[1])[j])
      );

      Rcpp::NumericVector mu_aux = auxParam[0];
      Rcpp::NumericVector nu_aux = auxParam[1];
      mu_aux.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      nu_aux.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      auxParam[0] = mu_aux;
      auxParam[1] = nu_aux;

      Rcpp::NumericVector lik = mixingDistribution->likelihood(data.row(i).t(), auxParam);
      probs[numberClusters + j] = (alpha / m) * lik[0];
    }

    // Normalize probabilities
    double probSum = 0.0;
    for (int j = 0; j < probs.size(); j++) {
      if (!std::isnan(probs[j]) && probs[j] >= 0) {
        probSum += probs[j];
      } else {
        probs[j] = 0.0;
      }
    }

    if (probSum <= 0) {
      // Uniform fallback
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

    // Update cluster assignment
    Rcpp::List updateResult = clusterLabelChange(i, newLabel, currentLabel, aux);

    // Update state from result
    clusterLabels = Rcpp::as<arma::uvec>(updateResult["clusterLabels"]);
    pointsPerCluster = Rcpp::as<arma::uvec>(updateResult["pointsPerCluster"]);
    clusterParameters = updateResult["clusterParameters"];
    numberClusters = updateResult["numberClusters"];
  }

  // Final validation
  arma::uword totalPoints = arma::sum(pointsPerCluster);
  if (totalPoints != n) {
    Rcpp::warning("Point count mismatch detected in C++ implementation");
  }
}

Rcpp::List NonConjugateBetaDP::clusterLabelChange(int i, int newLabel, int currentLabel,
                                                  const Rcpp::List& aux) {
  if (newLabel == currentLabel) {
    // Restore the point count since we temporarily removed it
    pointsPerCluster[currentLabel]++;
    return Rcpp::List::create(
      Rcpp::Named("clusterLabels") = clusterLabels,
      Rcpp::Named("pointsPerCluster") = pointsPerCluster,
      Rcpp::Named("clusterParameters") = clusterParameters,
      Rcpp::Named("numberClusters") = numberClusters
    );
  }

  // Extract current parameters
  Rcpp::NumericVector mu_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0]));
  Rcpp::NumericVector nu_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

  // Note: pointsPerCluster[currentLabel] has already been decremented in clusterComponentUpdate

  // Assign to new cluster
  if (newLabel < numberClusters) {
    // Existing cluster
    pointsPerCluster[newLabel]++;
    clusterLabels[i] = newLabel;

    // If old cluster is now empty, remove it
    if (pointsPerCluster[currentLabel] == 0 && currentLabel != newLabel) {
      numberClusters--;

      // Create new vectors without the empty cluster
      arma::uvec new_pointsPerCluster(numberClusters);
      int idx = 0;
      for (int j = 0; j < pointsPerCluster.n_elem; j++) {
        if (j != currentLabel) {
          new_pointsPerCluster[idx++] = pointsPerCluster[j];
        }
      }
      pointsPerCluster = new_pointsPerCluster;

      // Remove from parameter vectors
      Rcpp::NumericVector new_mu_vec;
      Rcpp::NumericVector new_nu_vec;
      for (int j = 0; j < mu_vec.size(); j++) {
        if (j != currentLabel) {
          new_mu_vec.push_back(mu_vec[j]);
          new_nu_vec.push_back(nu_vec[j]);
        }
      }
      mu_vec = new_mu_vec;
      nu_vec = new_nu_vec;

      // Update labels for points in clusters after the removed one
      for (arma::uword j = 0; j < clusterLabels.n_elem; j++) {
        if (clusterLabels[j] > (unsigned int)currentLabel) {
          clusterLabels[j]--;
        }
        // Also update the current point's label if needed
        if (j == i && newLabel > currentLabel) {
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
      Rcpp::NumericVector aux_nu = aux[1];

      mu_vec[currentLabel] = aux_mu[auxIndex];
      nu_vec[currentLabel] = aux_nu[auxIndex];

      pointsPerCluster[currentLabel] = 1;
      clusterLabels[i] = currentLabel;
    } else {
      // Create new cluster
      Rcpp::NumericVector aux_mu = aux[0];
      Rcpp::NumericVector aux_nu = aux[1];

      mu_vec.push_back(aux_mu[auxIndex]);
      nu_vec.push_back(aux_nu[auxIndex]);

      clusterLabels[i] = numberClusters;

      // Resize pointsPerCluster correctly
      arma::uvec new_pointsPerCluster(numberClusters + 1);
      for (int j = 0; j < numberClusters; j++) {
        new_pointsPerCluster[j] = pointsPerCluster[j];
      }
      new_pointsPerCluster[numberClusters] = 1;
      pointsPerCluster = new_pointsPerCluster;

      numberClusters++;
    }
  }

  // Update cluster parameters
  clusterParameters[0] = mu_vec;
  clusterParameters[1] = nu_vec;

  return Rcpp::List::create(
    Rcpp::Named("clusterLabels") = clusterLabels,
    Rcpp::Named("pointsPerCluster") = pointsPerCluster,
    Rcpp::Named("clusterParameters") = clusterParameters,
    Rcpp::Named("numberClusters") = numberClusters
  );
}

} // namespace dp
