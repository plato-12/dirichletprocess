// src/BetaDP.cpp
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"

namespace dp {

void NonConjugateBetaDP::clusterComponentUpdate() {
  // Implementation of Algorithm 8 from Neal (2000) for non-conjugate case
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
      Rcpp::NumericVector nu_vec = Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]);

      Rcpp::NumericVector mu_aux = aux[0];
      Rcpp::NumericVector nu_aux = aux[1];

      // Create new arrays including current cluster params
      Rcpp::NumericVector mu_combined(m);
      Rcpp::NumericVector nu_combined(m);
      mu_combined.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);
      nu_combined.attr("dim") = Rcpp::IntegerVector::create(1, 1, m);

      mu_combined[0] = mu_vec[currentLabel];
      nu_combined[0] = nu_vec[currentLabel];

      for (int j = 0; j < m-1; j++) {
        mu_combined[j+1] = mu_aux[j];
        nu_combined[j+1] = nu_aux[j];
      }

      aux = Rcpp::List::create(mu_combined, nu_combined);
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
        Rcpp::List clusterParam = Rcpp::List::create(
          Rcpp::NumericVector::create(Rcpp::as<Rcpp::NumericVector>(clusterParameters[0])[j]),
          Rcpp::NumericVector::create(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1])[j])
        );

        // Add dimensions
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

      // Add dimensions
      Rcpp::NumericVector mu_aux = auxParam[0];
      Rcpp::NumericVector nu_aux = auxParam[1];
      mu_aux.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      nu_aux.attr("dim") = Rcpp::IntegerVector::create(1, 1, 1);
      auxParam[0] = mu_aux;
      auxParam[1] = nu_aux;

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

    // Restore point count before update
    pointsPerCluster[currentLabel]++;

    // Update cluster assignment
    Rcpp::List updateResult = clusterLabelChange(i, newLabel, currentLabel, aux);

    // Update state from result
    clusterLabels = Rcpp::as<arma::uvec>(updateResult["clusterLabels"]);
    pointsPerCluster = Rcpp::as<arma::uvec>(updateResult["pointsPerCluster"]);
    clusterParameters = updateResult["clusterParameters"];
    numberClusters = updateResult["numberClusters"];
  }
}

Rcpp::List NonConjugateBetaDP::clusterLabelChange(int i, int newLabel, int currentLabel,
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
  Rcpp::NumericVector nu_vec = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(clusterParameters[1]));

  // 1. Remove point from old cluster
  pointsPerCluster[currentLabel]--;

  // 2. Assign to new cluster
  if (newLabel < numberClusters) {
    // Existing cluster
    pointsPerCluster[newLabel]++;
    clusterLabels[i] = newLabel;

    // If old cluster is now empty, remove it
    if (pointsPerCluster[currentLabel] == 0) {
      numberClusters--;
      pointsPerCluster.shed_row(currentLabel);

      mu_vec.erase(currentLabel);
      nu_vec.erase(currentLabel);

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
      pointsPerCluster.resize(numberClusters + 1);
      pointsPerCluster[numberClusters] = 1;
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
