// src/HierarchicalBetaDP.cpp
#include "../inst/include/HierarchicalDP.h"
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <RcppArmadillo.h>

namespace dp {

// HierarchicalBetaDP implementation
HierarchicalBetaDP::HierarchicalBetaDP() {
  // Constructor
}

HierarchicalBetaDP::~HierarchicalBetaDP() {
  // Clean up individual DPs
  for (auto& dp : indDP) {
    if (dp) {
      delete dp;
      dp = nullptr;
    }
  }
}

HierarchicalBetaDP* HierarchicalBetaDP::fromR(const Rcpp::List& rObj) {
  HierarchicalBetaDP* hdp = new HierarchicalBetaDP();

  try {
    if (rObj.containsElementNamed("indDP")) {
      Rcpp::List indDP_list = rObj["indDP"];

      for (int i = 0; i < indDP_list.size(); i++) {
        Rcpp::List dp_obj = indDP_list[i];

        // Create NonConjugateBetaDP from R object
        NonConjugateBetaDP* betaDP = new NonConjugateBetaDP();

        // Set common DP properties with validation
        if (dp_obj.containsElementNamed("data")) {
          betaDP->data = Rcpp::as<arma::mat>(dp_obj["data"]);
          betaDP->n = betaDP->data.n_rows;
        } else {
          delete betaDP;
          throw Rcpp::exception("Missing 'data' in DP object");
        }

        betaDP->alpha = dp_obj.containsElementNamed("alpha") ?
        Rcpp::as<double>(dp_obj["alpha"]) : 1.0;

        betaDP->alphaPriorParameters = dp_obj.containsElementNamed("alphaPriorParameters") ?
        Rcpp::as<Rcpp::NumericVector>(dp_obj["alphaPriorParameters"]) : Rcpp::NumericVector::create(1.0, 1.0);

        betaDP->mhDraws = dp_obj.containsElementNamed("mhDraws") ?
        Rcpp::as<int>(dp_obj["mhDraws"]) : 250;

        // Set cluster information with bounds checking
        if (dp_obj.containsElementNamed("clusterLabels")) {
          SEXP labels_sexp = dp_obj["clusterLabels"];
          if (!Rf_isNull(labels_sexp)) {
            arma::uvec labels = Rcpp::as<arma::uvec>(labels_sexp);
            if (labels.size() > 0) {
              // Labels should already be 0-indexed from R wrapper
              // Just check they're not negative
              if (labels.min() < 0) {
                delete betaDP;
                throw Rcpp::exception("Invalid cluster labels (must be >= 0 after conversion)");
              }
              betaDP->clusterLabels = labels;  // Already 0-indexed
            } else {
              // Empty labels vector
              betaDP->clusterLabels = arma::uvec();
            }
          } else {
            // NULL labels - initialize as empty
            betaDP->clusterLabels = arma::uvec();
          }
        } else {
          // No cluster labels - initialize as empty
          betaDP->clusterLabels = arma::uvec();
        }

        if (dp_obj.containsElementNamed("pointsPerCluster")) {
          betaDP->pointsPerCluster = Rcpp::as<arma::uvec>(dp_obj["pointsPerCluster"]);
        }

        betaDP->numberClusters = dp_obj.containsElementNamed("numberClusters") ?
        Rcpp::as<int>(dp_obj["numberClusters"]) : 1;

        if (dp_obj.containsElementNamed("clusterParameters")) {
          betaDP->clusterParameters = Rcpp::as<Rcpp::List>(dp_obj["clusterParameters"]);
        }

        betaDP->m = dp_obj.containsElementNamed("m") ?
        Rcpp::as<int>(dp_obj["m"]) : 3;

        // Create mixing distribution with validation
        if (dp_obj.containsElementNamed("mixingDistribution")) {
          Rcpp::List mixDist_obj = Rcpp::as<Rcpp::List>(dp_obj["mixingDistribution"]);

          if (mixDist_obj.containsElementNamed("priorParameters")) {
            betaDP->mixingDistribution = new BetaMixingDistribution(
              Rcpp::as<Rcpp::NumericVector>(mixDist_obj["priorParameters"]));

            betaDP->mixingDistribution->maxT = mixDist_obj.containsElementNamed("maxT") ?
            Rcpp::as<double>(mixDist_obj["maxT"]) : 1.0;

            if (mixDist_obj.containsElementNamed("mhStepSize")) {
              betaDP->mixingDistribution->mhStepSize = Rcpp::as<Rcpp::NumericVector>(mixDist_obj["mhStepSize"]);
            }
          } else {
            delete betaDP;
            throw Rcpp::exception("Missing prior parameters in mixing distribution");
          }
        } else {
          delete betaDP;
          throw Rcpp::exception("Missing mixing distribution");
        }

        hdp->indDP.push_back(betaDP);
      }
    }

    // Copy global parameters with validation
    if (rObj.containsElementNamed("globalParameters")) {
      hdp->globalParameters = Rcpp::as<Rcpp::List>(rObj["globalParameters"]);
    }

    if (rObj.containsElementNamed("globalStick")) {
      hdp->globalStick = Rcpp::as<arma::vec>(rObj["globalStick"]);
    }

    if (rObj.containsElementNamed("gamma")) {
      hdp->gamma = Rcpp::as<double>(rObj["gamma"]);
      if (hdp->gamma <= 0) {
        throw Rcpp::exception("gamma must be positive");
      }
    }

    if (rObj.containsElementNamed("gammaPriors")) {
      hdp->gammaPriors = Rcpp::as<Rcpp::NumericVector>(rObj["gammaPriors"]);
    }

    return hdp;

  } catch (std::exception& e) {
    // Clean up on error
    delete hdp;
    Rcpp::Rcerr << "Exception in HierarchicalBetaDP::fromR: " << e.what() << std::endl;
    throw;
  }
}

void HierarchicalBetaDP::fit(int iterations, bool updatePrior, bool progressBar) {
  if (progressBar) {
    Rcpp::Rcout << "Starting Hierarchical Beta DP fitting..." << std::endl;
  }

  // Store chain values by resizing the gammaChain member
  this->gammaChain = Rcpp::NumericVector(iterations);

  for (int iter = 0; iter < iterations; iter++) {
    // Update components
    clusterComponentUpdate();
    updateAlpha();
    globalParameterUpdate();
    updateG0();
    updateGamma();

    // Store gamma value
    this->gammaChain[iter] = gamma;

    // Update prior if requested
    if (updatePrior && indDP.size() > 0) {
      // Get all cluster parameters
      int total_clusters = 0;
      for (auto& dp : indDP) {
        NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(dp);
        if (betaDP) {
          total_clusters += betaDP->numberClusters;
        }
      }

      if (total_clusters > 0) {
        // Collect all nu parameters
        Rcpp::NumericVector all_nu;
        for (auto& dp : indDP) {
          NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(dp);
          if (betaDP) {
            Rcpp::List cluster_params_loop = betaDP->clusterParameters;
            if (cluster_params_loop.size() > 1) {
              Rcpp::NumericVector nu_params = Rcpp::as<Rcpp::NumericVector>(cluster_params_loop[1]);
              for (int i = 0; i < nu_params.size(); i++) {
                all_nu.push_back(nu_params[i]);
              }
            }
          }
        }

        // Update prior parameters using the first DP's mixing distribution
        if (indDP.empty() || !indDP[0]) {
          Rcpp::Rcerr << "indDP is empty or first element is null in fit()." << std::endl;
        } else {
          NonConjugateBetaDP* firstDP = dynamic_cast<NonConjugateBetaDP*>(indDP[0]);
          if (firstDP && firstDP->mixingDistribution) {
            Rcpp::List priorUpdateParams = Rcpp::List::create(
              Rcpp::Named("mu") = Rcpp::NumericVector(),
              Rcpp::Named("nu") = all_nu
            );
            firstDP->mixingDistribution->updatePriorParameters(priorUpdateParams, total_clusters);

            // Propagate updated prior to all DPs
            Rcpp::NumericVector newPrior = Rcpp::as<Rcpp::NumericVector>(
              firstDP->mixingDistribution->priorParameters);

            for (auto& dp_loop : indDP) {
              NonConjugateBetaDP* betaDP_loop = dynamic_cast<NonConjugateBetaDP*>(dp_loop);
              if (betaDP_loop && betaDP_loop->mixingDistribution) {
                betaDP_loop->mixingDistribution->priorParameters = newPrior;
              }
            }
          }
        }
      }
    }

    if (progressBar && (iterations == 0 || (iter + 1) % (iterations / 10) == 0 || iter == iterations - 1)) {
      Rcpp::Rcout << "Iteration " << iter + 1 << "/" << iterations << std::endl;
    }
  }

  if (progressBar) {
    Rcpp::Rcout << "Hierarchical Beta DP fitting complete." << std::endl;
  }
}

} // namespace dp
