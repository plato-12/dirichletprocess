// inst/include/mcmc_runner.h
#ifndef MCMC_RUNNER_H
#define MCMC_RUNNER_H

#include <RcppArmadillo.h>
#include <memory>
#include <vector>

namespace dirichletprocess {

// Forward declarations
class MixingDistribution;
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

  // Storage for results
  std::vector<arma::vec> alpha_samples;
  std::vector<std::vector<int>> cluster_samples;
  std::vector<std::vector<arma::vec>> theta_samples;

public:
  MCMCRunner(const arma::mat& data,
             const Rcpp::List& mixing_dist_params,
             const Rcpp::List& mcmc_params);

  // Main MCMC loop
  Rcpp::List run();

private:
  // MCMC steps
  void update_cluster_assignments();
  void update_cluster_parameters();
  void update_concentration();
  void store_iteration(int iter); // FIXED: Correct function name
};

// State container for DP
class DPState {
public:
  std::vector<int> cluster_labels;  // Use std::vector instead of arma::vec
  std::vector<arma::vec> cluster_params;
  arma::vec cluster_sizes;
  double alpha;
  int n_clusters;

  DPState(int n_obs, double initial_alpha)
    : cluster_labels(n_obs), alpha(initial_alpha), n_clusters(0) {
    cluster_sizes.set_size(0);
  }

  void update_cluster_counts() {
    // Ensure consistency
    n_clusters = cluster_params.size();
    if (n_clusters != static_cast<int>(cluster_sizes.n_elem)) {
      Rcpp::stop("Inconsistent cluster state");
    }
  }
};

} // namespace dirichletprocess

#endif
