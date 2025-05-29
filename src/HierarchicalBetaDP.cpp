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

    Rcpp::List clusterParamsList = betaDP->clusterParameters;
    if (clusterParamsList.size() == 0) continue;

    Rcpp::NumericVector mu_params = clusterParamsList[0];
    Rcpp::List globalParamsList = globalParameters;
    if (globalParamsList.size() == 0) continue;

    Rcpp::NumericVector mu_global = globalParamsList[0];

    // Bounds checking
    int n_clusters = std::min(betaDP->numberClusters, static_cast<int>(mu_params.size()));
    int n_global = mu_global.size();

    for (int j = 0; j < n_clusters; j++) {
      for (int k = 0; k < n_global; k++) {
        if (std::abs(mu_params[j] - mu_global[k]) < 1e-10) {
          all_global_labels.push_back(k);
          break;
        }
      }
    }
  }

  if (all_global_labels.empty()) return;

  std::sort(all_global_labels.begin(), all_global_labels.end());
  all_global_labels.erase(std::unique(all_global_labels.begin(), all_global_labels.end()),
                          all_global_labels.end());

  // Process each global parameter
  for (int global_idx : all_global_labels) {
    arma::mat combined_data;
    int total_points = 0;

    // Collect data points for this global parameter
    for (size_t dp_idx = 0; dp_idx < indDP.size(); dp_idx++) {
      NonConjugateBetaDP* betaDP = dynamic_cast<NonConjugateBetaDP*>(indDP[dp_idx]);
      if (!betaDP) continue;

      Rcpp::List clusterParamsList = betaDP->clusterParameters;
      if (clusterParamsList.size() == 0) continue;

      Rcpp::NumericVector mu_params = clusterParamsList[0];
      Rcpp::List globalParamsList = globalParameters;
      Rcpp::NumericVector mu_global = globalParamsList[0];

      // Bounds checking
      if (global_idx >= mu_global.size()) continue;

      int n_clusters = std::min(betaDP->numberClusters, static_cast<int>(mu_params.size()));

      for (int j = 0; j < n_clusters; j++) {
        if (std::abs(mu_params[j] - mu_global[global_idx]) < 1e-10) {
          // Check cluster labels bounds
          arma::uvec cluster_indices = arma::find(betaDP->clusterLabels == j);
          if (cluster_indices.n_elem > 0) {
            if (total_points == 0) {
              combined_data = betaDP->data.rows(cluster_indices);
            } else {
              combined_data = arma::join_vert(combined_data,
                                              betaDP->data.rows(cluster_indices));
            }
            total_points += cluster_indices.n_elem;
          }
        }
      }
    }

    if (total_points > 0 && indDP.size() > 0) {
      NonConjugateBetaDP* firstDP = dynamic_cast<NonConjugateBetaDP*>(indDP[0]);
      if (firstDP && firstDP->mixingDistribution) {
        BetaMixingDistribution* mixDist =
          dynamic_cast<BetaMixingDistribution*>(firstDP->getMixingDistribution());

        if (mixDist) {
          Rcpp::List new_params = mixDist->posteriorDraw(combined_data, 1);

          if (new_params.size() >= 2) {
            Rcpp::NumericVector new_mu = new_params[0];
            Rcpp::NumericVector new_nu = new_params[1];

            // Store old value before updating
            Rcpp::List globalParamsList = globalParameters;
            Rcpp::NumericVector mu_global_old = Rcpp::clone(
              Rcpp::as<Rcpp::NumericVector>(globalParamsList[0]));

            if (global_idx < mu_global_old.size()) {
              double old_mu_val = mu_global_old[global_idx];

              // Update global parameters
              Rcpp::NumericVector mu_global_new =
                Rcpp::as<Rcpp::NumericVector>(globalParamsList[0]);
              Rcpp::NumericVector nu_global_new =
                Rcpp::as<Rcpp::NumericVector>(globalParamsList[1]);

              if (new_mu.size() > 0 && new_nu.size() > 0) {
                mu_global_new[global_idx] = new_mu[0];
                nu_global_new[global_idx] = new_nu[0];
                globalParameters[0] = mu_global_new;
                globalParameters[1] = nu_global_new;

                // Update individual DP parameters
                for (size_t dp_idx_loop = 0; dp_idx_loop < indDP.size(); dp_idx_loop++) { // Renamed loop variable
                  NonConjugateBetaDP* betaDP_loop = // Renamed loop variable
                    dynamic_cast<NonConjugateBetaDP*>(indDP[dp_idx_loop]);
                  if (!betaDP_loop) continue;

                  Rcpp::List dpClusterParams = betaDP_loop->clusterParameters;
                  if (dpClusterParams.size() >= 2) {
                    Rcpp::NumericVector mu_params_dp = dpClusterParams[0];
                    Rcpp::NumericVector nu_params_dp = dpClusterParams[1];

                    for (int j = 0; j < mu_params_dp.size() && j < betaDP_loop->numberClusters; j++) {
                      if (std::abs(mu_params_dp[j] - old_mu_val) < 1e-10) {
                        mu_params_dp[j] = new_mu[0];
                        nu_params_dp[j] = new_nu[0];
                      }
                    }

                    betaDP_loop->clusterParameters[0] = mu_params_dp;
                    betaDP_loop->clusterParameters[1] = nu_params_dp;
                  }
                }
              }
            }
          }
        }
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
  // Ensure indDP is not empty before accessing its elements
  if (indDP.empty() || !indDP[0]) return; // Basic safety check
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

    // Ensure clusterParameters is not empty and has at least one element
    if (betaDP->clusterParameters.size() == 0) continue;
    Rcpp::NumericVector mu_params = betaDP->clusterParameters[0];

    // Ensure globalParameters is not empty and has at least one element
    if (globalParameters.size() == 0) continue;
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

  // Ensure gammaPriors has at least two elements
  if (gammaPriors.size() < 2) {
    // Handle error or set default, e.g., throw exception or return
    Rcpp::Rcerr << "gammaPriors does not have enough elements." << std::endl;
    return;
  }
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

  // Convert individual DPs with deep copy
  Rcpp::List indDP_list;
  for (auto& dp : indDP) {
    if (dp) {
      indDP_list.push_back(dp->toR());
    }
  }

  result["indDP"] = indDP_list;

  // Deep copy all global parameters
  result["globalParameters"] = Rcpp::clone(globalParameters);
  result["globalStick"] = Rcpp::clone(Rcpp::wrap(globalStick));
  result["gamma"] = gamma;
  result["gammaPriors"] = Rcpp::clone(gammaPriors);

  if (this->gammaChain.size() > 0) {
    result["gammaValues"] = Rcpp::NumericVector(this->gammaChain);
  }

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
          arma::uvec labels = Rcpp::as<arma::uvec>(dp_obj["clusterLabels"]);
          // Ensure labels are valid (>= 1 in R, >= 0 in C++ after conversion)
          if (labels.size() >0 && labels.min() < 1) { // Added labels.size() > 0 check
            delete betaDP;
            throw Rcpp::exception("Invalid cluster labels (must be >= 1)");
          }
          if (labels.size() > 0) { // Process only if labels exist
            betaDP->clusterLabels = labels - 1; // Convert to 0-indexed
          }
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
      // THE FIX IS APPLIED HERE:
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
        if (indDP.empty() || !indDP[0]) { // Safety check
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

            for (auto& dp_loop : indDP) { // Renamed loop variable
              NonConjugateBetaDP* betaDP_loop = dynamic_cast<NonConjugateBetaDP*>(dp_loop);
              if (betaDP_loop && betaDP_loop->mixingDistribution) {
                betaDP_loop->mixingDistribution->priorParameters = newPrior;
              }
            }
          }
        }
      }
    }

    if (progressBar && (iterations == 0 || (iter + 1) % (iterations / 10) == 0 || iter == iterations - 1)) { // Avoid division by zero
      Rcpp::Rcout << "Iteration " << iter + 1 << "/" << iterations << std::endl;
    }
  }

  if (progressBar) {
    Rcpp::Rcout << "Hierarchical Beta DP fitting complete." << std::endl;
  }
}

} // namespace dp
