// src/mcmc_runner.cpp
#include "../inst/include/mcmc_runner.h"
#include "../inst/include/mixing_distribution.h"
#include "../inst/include/utilities.h"
#include <algorithm>
#include <set>

namespace dirichletprocess {

MCMCRunner::MCMCRunner(const arma::mat& data,
                       const Rcpp::List& mixing_dist_params,
                       const Rcpp::List& mcmc_params)
  : data(data) {

  // Extract MCMC parameters
  n_iter = Rcpp::as<int>(mcmc_params["n_iter"]);
  n_burn = Rcpp::as<int>(mcmc_params["n_burn"]);
  thin = Rcpp::as<int>(mcmc_params["thin"]);
  update_concentration_flag = Rcpp::as<bool>(mcmc_params["update_concentration"]);

  // Set m_auxiliary (default to 3 for Algorithm 8)
  if (mcmc_params.containsElementNamed("m_auxiliary")) {
    m_auxiliary = Rcpp::as<int>(mcmc_params["m_auxiliary"]);
  } else {
    m_auxiliary = 3;
  }

  // Validate inputs
  if (data.n_rows < 2) {
    Rcpp::warning("Data has fewer than 2 observations. Results may be unreliable.");
  }

  // Create mixing distribution
  std::string dist_type = Rcpp::as<std::string>(mixing_dist_params["type"]);
  mixing_dist = MixingDistribution::create(dist_type, mixing_dist_params);

  // Initialize state
  double initial_alpha = Rcpp::as<double>(mcmc_params["alpha"]);
  if (initial_alpha <= 0) {
    Rcpp::stop("alpha must be positive");
  }

  state.reset(new DPState(data.n_rows, initial_alpha));

  // Pre-allocate storage
  alpha_samples.reserve(n_iter);
  cluster_samples.reserve(n_iter);
  theta_samples.reserve(n_iter);
}

Rcpp::List MCMCRunner::run() {
  // Initialize with one cluster containing all data
  std::fill(state->cluster_labels.begin(), state->cluster_labels.end(), 0);
  state->n_clusters = 1;
  state->cluster_sizes.set_size(1);
  state->cluster_sizes[0] = data.n_rows;

  // Initialize cluster parameters
  state->cluster_params.resize(1);
  state->cluster_params[0] = mixing_dist->prior_draw();

  // MCMC loop
  for (int iter = 0; iter < n_iter; ++iter) {
    // Update cluster assignments using Algorithm 8
    update_cluster_assignments_algorithm8();

    // Update cluster parameters
    update_cluster_parameters();

    // Update concentration parameter
    if (update_concentration_flag) {
      update_concentration();
    }

    // Store current iteration
    store_iteration(iter);
  }

  // Convert stored samples to proper format for R
  int n_stored = 0;
  std::vector<arma::vec> alpha_chain;
  std::vector<std::vector<int>> labels_chain;
  std::vector<std::vector<arma::vec>> theta_chain;
  std::vector<int> n_clusters_chain;

  for (int iter = n_burn; iter < n_iter; iter += thin) {
    alpha_chain.push_back(alpha_samples[iter]);
    labels_chain.push_back(cluster_samples[iter]);
    theta_chain.push_back(theta_samples[iter]);

    std::set<int> unique_labels(cluster_samples[iter].begin(),
                                cluster_samples[iter].end());
    n_clusters_chain.push_back(unique_labels.size());
    n_stored++;
  }

  // Convert to matrices for R
  arma::mat labels_matrix(n_stored, data.n_rows);
  arma::vec alpha_vector(n_stored);

  for (int i = 0; i < n_stored; ++i) {
    alpha_vector[i] = alpha_chain[i][0];
    for (arma::uword j = 0; j < data.n_rows; j++) {
      labels_matrix(i, j) = labels_chain[i][j] + 1; // Convert to 1-indexed
    }
  }

  // Get final state
  int final_iter = n_iter - 1;
  arma::vec final_labels = arma::conv_to<arma::vec>::from(cluster_samples[final_iter]) + 1;

  // Convert final theta to List
  Rcpp::List final_theta_list(theta_samples[final_iter].size());
  for (size_t i = 0; i < theta_samples[final_iter].size(); ++i) {
    final_theta_list[i] = theta_samples[final_iter][i];
  }

  // Convert theta_chain to List
  Rcpp::List theta_chain_list(n_stored);
  for (int i = 0; i < n_stored; ++i) {
    Rcpp::List iter_params(theta_chain[i].size());
    for (size_t j = 0; j < theta_chain[i].size(); ++j) {
      iter_params[j] = theta_chain[i][j];
    }
    theta_chain_list[i] = iter_params;
  }

  return Rcpp::List::create(
    Rcpp::Named("cluster_labels") = labels_matrix,
    Rcpp::Named("alpha") = alpha_vector,
    Rcpp::Named("theta") = final_theta_list,
    Rcpp::Named("n_clusters") = n_clusters_chain,
    Rcpp::Named("final_labels") = final_labels,
    Rcpp::Named("final_n_clusters") = state->n_clusters,
    Rcpp::Named("alpha_chain") = alpha_vector,
    Rcpp::Named("labels_chain") = labels_matrix,
    Rcpp::Named("theta_chain") = theta_chain_list
  );
}

void MCMCRunner::update_cluster_assignments_algorithm8() {
  int n = data.n_rows;

  for (int i = 0; i < n; ++i) {
    arma::vec obs = data.row(i).t();
    int current_cluster = state->cluster_labels[i];

    // Remove observation from current cluster
    state->cluster_sizes[current_cluster]--;

    // Collect existing cluster parameters
    std::vector<arma::vec> existing_params;
    std::vector<int> existing_indices;

    if (state->cluster_sizes[current_cluster] == 0) {
      // Current cluster is now empty - include it as one of the auxiliary parameters
      existing_params.push_back(state->cluster_params[current_cluster]);
      existing_indices.push_back(current_cluster);

      // Draw m-1 additional auxiliary parameters
      for (int j = 0; j < m_auxiliary - 1; ++j) {
        existing_params.push_back(mixing_dist->prior_draw());
        existing_indices.push_back(-1); // Mark as auxiliary
      }
    } else {
      // Current cluster still has points
      for (int k = 0; k < state->n_clusters; ++k) {
        if (state->cluster_sizes[k] > 0) {
          existing_params.push_back(state->cluster_params[k]);
          existing_indices.push_back(k);
        }
      }

      // Draw m auxiliary parameters
      for (int j = 0; j < m_auxiliary; ++j) {
        existing_params.push_back(mixing_dist->prior_draw());
        existing_indices.push_back(-1); // Mark as auxiliary
      }
    }

    // Calculate probabilities
    arma::vec probs(existing_params.size());

    for (size_t j = 0; j < existing_params.size(); ++j) {
      double log_lik = mixing_dist->log_likelihood(obs, existing_params[j]);

      if (existing_indices[j] >= 0) {
        // Existing cluster
        probs[j] = state->cluster_sizes[existing_indices[j]] * std::exp(log_lik);
      } else {
        // Auxiliary cluster
        probs[j] = (state->alpha / m_auxiliary) * std::exp(log_lik);
      }
    }

    // Normalize probabilities
    double prob_sum = arma::sum(probs);
    if (prob_sum > 0) {
      probs = probs / prob_sum;
    } else {
      probs.fill(1.0 / probs.n_elem);
    }

    // Sample new cluster
    int chosen_idx = sample_categorical(probs);

    if (existing_indices[chosen_idx] >= 0) {
      // Assign to existing cluster
      state->cluster_labels[i] = existing_indices[chosen_idx];
      state->cluster_sizes[existing_indices[chosen_idx]]++;
    } else {
      // Create new cluster with the chosen auxiliary parameters
      if (state->cluster_sizes[current_cluster] == 0) {
        // Reuse the empty cluster slot
        state->cluster_labels[i] = current_cluster;
        state->cluster_sizes[current_cluster] = 1;
        state->cluster_params[current_cluster] = existing_params[chosen_idx];
      } else {
        // Add new cluster
        state->cluster_labels[i] = state->n_clusters;
        state->cluster_params.push_back(existing_params[chosen_idx]);
        state->cluster_sizes.resize(state->n_clusters + 1);
        state->cluster_sizes[state->n_clusters] = 1;
        state->n_clusters++;
      }
    }
  }

  // Clean up empty clusters and relabel
  std::vector<int> new_labels(state->n_clusters, -1);
  int new_idx = 0;

  for (int k = 0; k < state->n_clusters; ++k) {
    if (state->cluster_sizes[k] > 0) {
      new_labels[k] = new_idx++;
    }
  }

  // Update cluster labels and parameters
  std::vector<arma::vec> new_params;
  for (int k = 0; k < state->n_clusters; ++k) {
    if (new_labels[k] >= 0) {
      new_params.push_back(state->cluster_params[k]);
    }
  }

  // Relabel data points
  for (int i = 0; i < n; ++i) {
    state->cluster_labels[i] = new_labels[state->cluster_labels[i]];
  }

  state->cluster_params = new_params;
  state->n_clusters = new_params.size();
  state->update_cluster_counts();
}

void MCMCRunner::update_cluster_parameters() {
  for (int k = 0; k < state->n_clusters; ++k) {
    // Get data points in cluster k
    std::vector<arma::uword> cluster_indices;

    for (arma::uword i = 0; i < data.n_rows; i++) {
      if (state->cluster_labels[i] == k) {
        cluster_indices.push_back(i);
      }
    }

    if (cluster_indices.size() > 0) {
      // Extract cluster data
      arma::mat cluster_data(cluster_indices.size(), data.n_cols);

      for (size_t idx = 0; idx < cluster_indices.size(); ++idx) {
        cluster_data.row(idx) = data.row(cluster_indices[idx]);
      }

      // Draw from posterior
      state->cluster_params[k] = mixing_dist->posterior_draw(cluster_data, state->cluster_params[k]);
    }
  }
}

void MCMCRunner::update_concentration() {
  // Auxiliary variable method for updating alpha
  double a = 1.0;  // Gamma prior shape
  double b = 1.0;  // Gamma prior rate

  double eta = R::rbeta(state->alpha + 1, data.n_rows);
  double pi_eta = (a + state->n_clusters - 1) /
    (a + state->n_clusters - 1 + data.n_rows * (b - std::log(eta)));

  if (R::runif(0, 1) < pi_eta) {
    state->alpha = R::rgamma(a + state->n_clusters, 1.0 / (b - std::log(eta)));
  } else {
    state->alpha = R::rgamma(a + state->n_clusters - 1, 1.0 / (b - std::log(eta)));
  }
}

void MCMCRunner::store_iteration(int iter) {
  alpha_samples.push_back(arma::vec{state->alpha});

  std::vector<int> labels(state->cluster_labels.begin(), state->cluster_labels.end());
  cluster_samples.push_back(labels);

  std::vector<arma::vec> params_copy;
  for (const auto& param : state->cluster_params) {
    params_copy.push_back(param);
  }
  theta_samples.push_back(params_copy);
}

} // namespace dirichletprocess
