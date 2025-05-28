// src/HierarchicalBetaDP.cpp
#include "../inst/include/HierarchicalDP.h"
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <RcppArmadillo.h>

namespace dp {

// HierarchicalDP base class implementation
HierarchicalDP::HierarchicalDP() : gamma(1.0) {
  // Constructor
}

HierarchicalDP::~HierarchicalDP() {
  // Destructor - cleanup is handled in derived classes
}

void HierarchicalDP::clusterComponentUpdate() {
  // Update each individual DP
  for (auto& dp : indDP) {
    if (dp) {
      dp->clusterComponentUpdate();
    }
  }
}

void HierarchicalDP::clusterParameterUpdate() {
  // Update parameters for each individual DP
  for (auto& dp : indDP) {
    if (dp) {
      dp->clusterParameterUpdate();
    }
  }
}

void HierarchicalDP::updateAlpha() {
  // Update alpha for each individual DP
  for (auto& dp : indDP) {
    if (dp) {
      dp->updateAlpha();
    }
  }
}

void HierarchicalDP::globalParameterUpdate() {
  // Get unique global labels across all DPs
  std::vector<int> all_global_labels;

  for (size_t i = 0; i < indDP.size(); i++) {
    NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(indDP[i]);
    if (!betaDP) continue;

    Rcpp::NumericVector mu_params = Rcpp::as<Rcpp::List>(betaDP->clusterParameters)[0];
    Rcpp::NumericVector mu_global = Rcpp::as<Rcpp::List>(globalParameters)[0];

    for (int j = 0; j < betaDP->numberClusters; j++) {
      for (int k = 0; k < mu_global.size(); k++) {
        if (std::abs(mu_params[j] - mu_global[k]) < 1e-10) {
          all_global_labels.push_back(k);
          break;
        }
      }
    }
  }

  std::sort(all_global_labels.begin(), all_global_labels.end());
  all_global_labels.erase(std::unique(all_global_labels.begin(), all_global_labels.end()), all_global_labels.end());

  for (int global_idx : all_global_labels) {
    arma::mat combined_data;
    int total_points = 0;

    for (size_t dp_idx = 0; dp_idx < indDP.size(); dp_idx++) {
      NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(indDP[dp_idx]);
      if (!betaDP) continue;

      Rcpp::NumericVector mu_params = Rcpp::as<Rcpp::List>(betaDP->clusterParameters)[0];
      Rcpp::NumericVector mu_global = Rcpp::as<Rcpp::List>(globalParameters)[0];

      for (int j = 0; j < betaDP->numberClusters; j++) {
        if (std::abs(mu_params[j] - mu_global[global_idx]) < 1e-10) {
          arma::uvec cluster_indices = arma::find(betaDP->clusterLabels == j);
          if (cluster_indices.n_elem > 0) {
            if (total_points == 0) {
              combined_data = betaDP->data.rows(cluster_indices);
            } else {
              combined_data = arma::join_vert(combined_data, betaDP->data.rows(cluster_indices));
            }
            total_points += cluster_indices.n_elem;
          }
        }
      }
    }

    if (total_points > 0) {
      BetaMixingDistribution* mixDist = dynamic_cast<BetaMixingDistribution*>(indDP[0]->getMixingDistribution());
      if (mixDist) {
        Rcpp::List new_params = mixDist->posteriorDraw(combined_data, 1);
        Rcpp::NumericVector new_mu = new_params[0];
        Rcpp::NumericVector new_nu = new_params[1];

        // *** START OF CRITICAL FIX ***
        // Store the old global parameter value before updating
        Rcpp::NumericVector mu_global_old = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(globalParameters[0]));
        double old_mu_val = mu_global_old[global_idx];

        // Update global parameters
        Rcpp::NumericVector mu_global_new = Rcpp::as<Rcpp::NumericVector>(globalParameters[0]);
        Rcpp::NumericVector nu_global_new = Rcpp::as<Rcpp::NumericVector>(globalParameters[1]);
        mu_global_new[global_idx] = new_mu[0];
        nu_global_new[global_idx] = new_nu[0];
        globalParameters[0] = mu_global_new;
        globalParameters[1] = nu_global_new;

        // Update individual DP parameters that were using the old global parameter
        for (size_t dp_idx = 0; dp_idx < indDP.size(); dp_idx++) {
          NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(indDP[dp_idx]);
          if (!betaDP) continue;

          Rcpp::NumericVector mu_params = betaDP->clusterParameters[0];
          Rcpp::NumericVector nu_params = betaDP->clusterParameters[1];

          for (int j = 0; j < betaDP->numberClusters; j++) {
            // Use the stored old value for comparison
            if (std::abs(mu_params[j] - old_mu_val) < 1e-10) {
              mu_params[j] = new_mu[0];
              nu_params[j] = new_nu[0];
            }
          }

          betaDP->clusterParameters[0] = mu_params;
          betaDP->clusterParameters[1] = nu_params;
        }
        // *** END OF CRITICAL FIX ***
      }
    }
  }
}

void HierarchicalDP::updateG0() {
  // Get global parameters and their frequencies
  std::map<int, int> global_param_counts;

  for (size_t i = 0; i < indDP.size(); i++) {
    NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(indDP[i]);
    if (!betaDP) continue;

    Rcpp::NumericVector mu_params = betaDP->clusterParameters[0];
    Rcpp::NumericVector mu_global = globalParameters[0];

    for (int j = 0; j < betaDP->numberClusters; j++) {
      for (int k = 0; k < mu_global.size(); k++) {
        if (std::abs(mu_params[j] - mu_global[k]) < 1e-10) {
          global_param_counts[k]++;
          break;
        }
      }
    }
  }

  int num_tables = global_param_counts.size();
  if (num_tables == 0) return;

  // Get frequencies
  Rcpp::NumericVector frequencies(num_tables);
  int idx = 0;
  for (auto& pair : global_param_counts) {
    frequencies[idx++] = pair.second;
  }

  // Draw from Dirichlet distribution
  Rcpp::NumericVector dirichlet_params = Rcpp::NumericVector::create();
  for (int i = 0; i < num_tables; i++) {
    dirichlet_params.push_back(frequencies[i]);
  }
  dirichlet_params.push_back(gamma);

  // Use R's rdirichlet through Rcpp
  Rcpp::Environment gtools("package:gtools");
  Rcpp::Function rdirichlet = gtools["rdirichlet"];
  Rcpp::NumericMatrix dirichlet_draw = rdirichlet(1, dirichlet_params);
  Rcpp::NumericVector weights = dirichlet_draw(0, Rcpp::_);

  // Update stick breaking weights
  int num_breaks = std::ceil(gamma + num_tables) * 20 + 5;
  globalStick.set_size(num_tables + num_breaks);

  // Existing table weights
  for (int i = 0; i < num_tables; i++) {
    globalStick[i] = weights[i];
  }

  // New table weights from stick breaking
  double remaining_weight = weights[num_tables];
  for (int i = 0; i < num_breaks; i++) {
    double beta = R::rbeta(1.0, gamma + num_tables);
    globalStick[num_tables + i] = beta * remaining_weight;
    remaining_weight *= (1.0 - beta);
  }

  // Draw new parameters for the additional breaks
  BetaMixingDistribution* mixDist = dynamic_cast<BetaMixingDistribution*>(
    indDP[0]->getMixingDistribution());

  if (mixDist) {
    Rcpp::List new_params = mixDist->priorDraw(num_breaks);

    // Expand global parameters
    Rcpp::NumericVector mu_global = globalParameters[0];
    Rcpp::NumericVector nu_global = globalParameters[1];
    Rcpp::NumericVector new_mu = new_params[0];
    Rcpp::NumericVector new_nu = new_params[1];

    for (int i = 0; i < num_breaks; i++) {
      mu_global.push_back(new_mu[i]);
      nu_global.push_back(new_nu[i]);
    }

    globalParameters[0] = mu_global;
    globalParameters[1] = nu_global;
  }
}

void HierarchicalDP::updateGamma() {
  // Count unique global parameters
  std::set<int> unique_params;

  for (size_t i = 0; i < indDP.size(); i++) {
    NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(indDP[i]);
    if (!betaDP) continue;

    Rcpp::NumericVector mu_params = betaDP->clusterParameters[0];
    Rcpp::NumericVector mu_global = globalParameters[0];

    for (int j = 0; j < betaDP->numberClusters; j++) {
      for (int k = 0; k < mu_global.size(); k++) {
        if (std::abs(mu_params[j] - mu_global[k]) < 1e-10) {
          unique_params.insert(k);
          break;
        }
      }
    }
  }

  int num_unique = unique_params.size();
  int num_tables = 0;

  for (auto& dp : indDP) {
    NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(dp);
    if (betaDP) {
      num_tables += betaDP->numberClusters;
    }
  }

  // Update gamma using auxiliary variable method
  double x = R::rbeta(gamma + 1.0, num_tables);
  double log_x = std::log(x);

  double pi1 = gammaPriors[0] + num_unique - 1.0;
  double pi2 = num_tables * (gammaPriors[1] - log_x);

  double pi_ratio = pi1 / (pi1 + pi2);
  if (!std::isfinite(pi_ratio)) {
    pi_ratio = 0.5;
  }

  double post_shape;
  if (R::runif(0, 1) < pi_ratio) {
    post_shape = gammaPriors[0] + num_unique;
  } else {
    post_shape = gammaPriors[0] + num_unique - 1.0;
  }

  double post_rate = gammaPriors[1] - log_x;
  if (post_rate <= 0) post_rate = 1e-6;

  gamma = R::rgamma(post_shape, 1.0 / post_rate);
  if (gamma <= 0) gamma = 1e-6;
}

Rcpp::List HierarchicalDP::toR() const {
  Rcpp::List result;

  // Convert individual DPs
  Rcpp::List indDP_list;
  for (auto& dp : indDP) {
    if (dp) {
      indDP_list.push_back(dp->toR());
    }
  }

  result["indDP"] = indDP_list;
  result["globalParameters"] = globalParameters;
  result["globalStick"] = Rcpp::wrap(globalStick);
  result["gamma"] = gamma;
  result["gammaPriors"] = gammaPriors;
  result["gammaValues"] = this->gammaChain;

  result.attr("class") = Rcpp::CharacterVector::create("list", "dirichletprocess", "hierarchical");

  return result;
}

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

  if (rObj.containsElementNamed("indDP")) {
    Rcpp::List indDP_list = rObj["indDP"];

    for (int i = 0; i < indDP_list.size(); i++) {
      Rcpp::List dp_obj = indDP_list[i];

      // Create NonConjugateBetaDP from R object
      NonConjugateBetaDP* betaDP = new NonConjugateBetaDP();

      // Set common DP properties
      betaDP->data = Rcpp::as<arma::mat>(dp_obj["data"]);
      betaDP->n = betaDP->data.n_rows;
      betaDP->alpha = Rcpp::as<double>(dp_obj["alpha"]);
      betaDP->alphaPriorParameters = dp_obj["alphaPriorParameters"];
      betaDP->mhDraws = dp_obj.containsElementNamed("mhDraws") ?
      Rcpp::as<int>(dp_obj["mhDraws"]) : 250;

      // Set cluster information
      betaDP->clusterLabels = Rcpp::as<arma::uvec>(dp_obj["clusterLabels"]) - 1; // Convert to 0-indexed
      betaDP->pointsPerCluster = Rcpp::as<arma::uvec>(dp_obj["pointsPerCluster"]);
      betaDP->numberClusters = Rcpp::as<int>(dp_obj["numberClusters"]);
      betaDP->clusterParameters = dp_obj["clusterParameters"];
      betaDP->m = dp_obj.containsElementNamed("m") ? Rcpp::as<int>(dp_obj["m"]) : 3;

      // Create mixing distribution
      Rcpp::List mixDist = dp_obj["mixingDistribution"];
      betaDP->mixingDistribution = new BetaMixingDistribution(
        Rcpp::as<Rcpp::NumericVector>(mixDist["priorParameters"]));
      betaDP->mixingDistribution->maxT = mixDist.containsElementNamed("maxT") ?
      Rcpp::as<double>(mixDist["maxT"]) : 1.0;

      if (mixDist.containsElementNamed("mhStepSize")) {
        betaDP->mixingDistribution->mhStepSize = mixDist["mhStepSize"];
      }

      hdp->indDP.push_back(betaDP);
    }
  }

  if (rObj.containsElementNamed("globalParameters")) {
    hdp->globalParameters = rObj["globalParameters"];
  }

  if (rObj.containsElementNamed("globalStick")) {
    hdp->globalStick = Rcpp::as<arma::vec>(rObj["globalStick"]);
  }

  if (rObj.containsElementNamed("gamma")) {
    hdp->gamma = Rcpp::as<double>(rObj["gamma"]);
  }

  if (rObj.containsElementNamed("gammaPriors")) {
    hdp->gammaPriors = Rcpp::as<Rcpp::NumericVector>(rObj["gammaPriors"]);
  }

  return hdp;
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
            Rcpp::NumericVector nu_params = betaDP->clusterParameters[1];
            for (int i = 0; i < nu_params.size(); i++) {
              all_nu.push_back(nu_params[i]);
            }
          }
        }

        // Update prior parameters using the first DP's mixing distribution
        NonConjugateBetaDP* firstDP = dynamic_cast<NonConjugateBetaDP*>(indDP[0]);
        if (firstDP && firstDP->mixingDistribution) {
          Rcpp::List clusterParams = Rcpp::List::create(
            Rcpp::Named("mu") = Rcpp::NumericVector(),
            Rcpp::Named("nu") = all_nu
          );
          firstDP->mixingDistribution->updatePriorParameters(clusterParams, total_clusters);

          // Propagate updated prior to all DPs
          Rcpp::NumericVector newPrior = Rcpp::as<Rcpp::NumericVector>(
            firstDP->mixingDistribution->priorParameters);

          for (auto& dp : indDP) {
            NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(dp);
            if (betaDP && betaDP->mixingDistribution) {
              betaDP->mixingDistribution->priorParameters = newPrior;
            }
          }
        }
      }
    }

    if (progressBar && ((iter + 1) % (iterations / 10) == 0 || iter == iterations - 1)) {
      Rcpp::Rcout << "Iteration " << iter + 1 << "/" << iterations << std::endl;
    }
  }

  if (progressBar) {
    Rcpp::Rcout << "Hierarchical Beta DP fitting complete." << std::endl;
  }
}

} // namespace dp
