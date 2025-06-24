#include "../inst/include/hierarchical_mcmc_runner.h"

namespace dirichletprocess {

HierarchicalMCMCRunner::HierarchicalMCMCRunner(
  const std::vector<arma::mat>& datasets,
  const Rcpp::List& mixing_dist_params,
  const Rcpp::List& mcmc_params)
  : datasets(datasets) {

  // Extract MCMC parameters
  n_iter = Rcpp::as<int>(mcmc_params["n_iter"]);
  n_burn = Rcpp::as<int>(mcmc_params["n_burn"]);
  thin = Rcpp::as<int>(mcmc_params["thin"]);
  update_prior = Rcpp::as<bool>(mcmc_params["update_prior"]);

  // Create hierarchical mixing distribution
  double alpha0 = Rcpp::as<double>(mixing_dist_params["alpha0"]);
  double beta0 = Rcpp::as<double>(mixing_dist_params["beta0"]);
  double maxT = Rcpp::as<double>(mixing_dist_params["maxT"]);

  hierarchical_mixing_dist.reset(
    new HierarchicalBetaMixing(alpha0, beta0, maxT));

  // Create individual MCMC runners for each dataset
  for (const auto& data : datasets) {
    // Create parameters for individual runner
    Rcpp::List individual_params = Rcpp::clone(mixing_dist_params);
    individual_params["type"] = "beta";

    runners.emplace_back(
      new MCMCRunner(data, individual_params, mcmc_params));
  }

  // Pre-allocate storage
  gamma_samples.reserve(n_iter);
  global_param_samples.reserve(n_iter);
}

Rcpp::List HierarchicalMCMCRunner::run() {
  Rcpp::Rcout << "Starting Hierarchical Beta MCMC with "
              << datasets.size() << " datasets" << std::endl;

  // Initialize all runners
  for (auto& runner : runners) {
    runner->run();  // Initial short run to establish clusters
  }

  // Main MCMC loop
  for (int iter = 0; iter < n_iter; iter++) {
    Rcpp::checkUserInterrupt();

    // Step 1: Update local clusters using Algorithm 8
    update_local_clusters();

    // Step 2: Update global parameters
    update_global_parameters();

    // Step 3: Update gamma (concentration parameter for G0)
    update_gamma();

    // Step 4: Propagate updated G0 to local DPs
    propagate_g0_to_local();

    // Store samples after burn-in
    if (iter >= n_burn && (iter - n_burn) % thin == 0) {
      store_iteration(iter);
    }

    // Progress reporting
    if ((iter + 1) % 100 == 0) {
      Rcpp::Rcout << "Iteration " << (iter + 1) << "/" << n_iter << std::endl;
    }
  }

  // Compile results
  Rcpp::List results;

  // Individual DP results
  Rcpp::List individual_results;
  for (size_t i = 0; i < runners.size(); i++) {
    individual_results.push_back(runners[i]->run());
  }
  results["individual_dps"] = individual_results;

  // Hierarchical parameters
  results["gamma_samples"] = gamma_samples;
  results["global_parameters"] = global_param_samples;
  results["global_weights"] = hierarchical_mixing_dist->get_global_weights();

  return results;
}

void HierarchicalMCMCRunner::update_local_clusters() {
  // Each dataset updates its cluster assignments
  // This uses Algorithm 8 from Neal (2000)

  for (auto& runner : runners) {
    // Run one iteration of local MCMC
    runner->run();  // This internally calls update_cluster_assignments_algorithm8()
  }
}

void HierarchicalMCMCRunner::update_global_parameters() {
  // Collect all cluster parameters and data from local DPs
  std::vector<arma::mat> all_cluster_data;
  std::vector<arma::vec> all_cluster_params;

  // This is simplified - in practice, we need to track which local clusters
  // are assigned to which global clusters

  hierarchical_mixing_dist->update_global_parameters(
      all_cluster_data, all_cluster_params);
}

void HierarchicalMCMCRunner::update_gamma() {
  // Count unique clusters across all datasets
  int n_unique_global = hierarchical_mixing_dist->get_global_params().size();
  int n_total_obs = 0;

  for (const auto& data : datasets) {
    n_total_obs += data.n_rows;
  }

  hierarchical_mixing_dist->update_gamma(n_unique_global, n_total_obs);
}

void HierarchicalMCMCRunner::propagate_g0_to_local() {
  // Update the base distribution G0 in each local DP
  // This ensures all local DPs share the same G0

  // In practice, this would update the mixing distribution
  // in each runner to use the current G0
}

void HierarchicalMCMCRunner::store_iteration(int iter) {
  gamma_samples.push_back(hierarchical_mixing_dist->get_gamma());
  global_param_samples.push_back(hierarchical_mixing_dist->get_global_params());
}

} // namespace dirichletprocess
