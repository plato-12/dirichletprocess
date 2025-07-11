// inst/include/mcmc_runner.h
#ifndef MCMC_RUNNER_H
#define MCMC_RUNNER_H

#include <RcppArmadillo.h>
#include <memory>
#include <vector>
#include <map>
#include "mixing_distribution_base.h"  // CHANGE: Include full header instead of forward declaration

namespace dirichletprocess {

// Forward declarations
// class MixingDistribution;  // REMOVE THIS LINE - not needed anymore
class DPState;

// Main MCMC runner class
class MCMCRunner {
private:
  // Data
  arma::mat data;

  // Model components
  std::unique_ptr<MixingDistribution> mixing_dist;

  // MCMC state
  std::unique_ptr<DPState> state;

  // Parameters
  int n_iter;
  int n_burn;
  int thin;
  bool update_concentration_flag;
  int m_auxiliary; // Number of auxiliary parameters for Algorithm 8

  // Alpha prior parameters
  double alpha_prior_shape;
  double alpha_prior_rate;

  // Storage for results
  std::vector<arma::vec> alpha_samples;
  std::vector<std::vector<int>> cluster_samples;
  std::vector<std::vector<arma::vec>> theta_samples;
  std::vector<double> likelihood_samples;  // For tracking likelihood

public:
  MCMCRunner(const arma::mat& data,
             const Rcpp::List& mixing_dist_params,
             const Rcpp::List& mcmc_params);

  // Main MCMC loop
  Rcpp::List run();

private:
  // MCMC steps
  void update_cluster_assignments_algorithm8(); // Algorithm 8 implementation
  void update_cluster_parameters();
  void update_concentration();
  void store_iteration(int iter);
  void cleanup_empty_clusters();

  // Helper function for categorical sampling
  int sample_categorical(const std::vector<double>& probs);
};

// State container for DP
class DPState {
public:
  std::vector<int> cluster_labels;
  std::vector<arma::vec> cluster_params;
  arma::vec cluster_sizes;
  double alpha;
  int n_clusters;

  DPState(int n_obs, double initial_alpha)
    : cluster_labels(n_obs), alpha(initial_alpha), n_clusters(0) {
    cluster_sizes.set_size(0);
  }

  void update_cluster_counts() {
    // Count actual clusters
    std::map<int, int> label_counts;
    for (int label : cluster_labels) {
      label_counts[label]++;
    }

    n_clusters = label_counts.size();

    // Update cluster sizes
    cluster_sizes.zeros(n_clusters);
    int idx = 0;
    for (const auto& pair : label_counts) {
      cluster_sizes[idx] = pair.second;
      idx++;
    }
  }
};

} // namespace dirichletprocess

#endif
