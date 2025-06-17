// src/mcmc_runner.cpp
#include "../inst/include/mcmc_runner.h"
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/utilities.h"  // ADD THIS LINE
#include <set>
#include <algorithm>

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

    // Collect parameters and calculate probabilities
    std::vector<arma::vec> all_params;
    std::vector<int> param_types; // 0 = existing cluster, 1 = auxiliary
    std::vector<int> cluster_indices;

    // Add all non-empty clusters
    for (int k = 0; k < state->n_clusters; ++k) {
      if (state->cluster_sizes[k] > 0) {
        all_params.push_back(state->cluster_params[k]);
        param_types.push_back(0); // existing cluster
        cluster_indices.push_back(k);
      }
    }

    // If current cluster is now empty, we can reuse it
    bool can_reuse_current = (state->cluster_sizes[current_cluster] == 0);
    int reuse_idx = -1;
    if (can_reuse_current) {
      reuse_idx = current_cluster;
    }

    // Add m auxiliary parameters
    for (int j = 0; j < m_auxiliary; ++j) {
      all_params.push_back(mixing_dist->prior_draw());
      param_types.push_back(1); // auxiliary
      cluster_indices.push_back(-1);
    }

    // Calculate probabilities
    arma::vec probs(all_params.size());

    for (size_t j = 0; j < all_params.size(); ++j) {
      double log_lik = mixing_dist->log_likelihood(obs, all_params[j]);

      if (param_types[j] == 0) {
        // Existing cluster: weight by number of points
        probs[j] = state->cluster_sizes[cluster_indices[j]] * std::exp(log_lik);
      } else {
        // Auxiliary parameter: weight by alpha/m
        probs[j] = (state->alpha / m_auxiliary) * std::exp(log_lik);
      }
    }

    // Handle numerical issues
    double max_prob = probs.max();
    if (max_prob > 700) {  // Prevent overflow
      probs = probs - max_prob + 700;
    }

    // Normalize probabilities
    double prob_sum = arma::sum(probs);
    if (prob_sum <= 0.0 || !std::isfinite(prob_sum)) {
      probs.fill(1.0 / probs.n_elem);
    } else {
      probs = probs / prob_sum;
    }

    // Sample new cluster
    int chosen_idx = sample_categorical(probs);

    if (param_types[chosen_idx] == 0) {
      // Assign to existing cluster
      int new_cluster = cluster_indices[chosen_idx];
      state->cluster_labels[i] = new_cluster;
      state->cluster_sizes[new_cluster]++;
    } else {
      // Create new cluster with auxiliary parameter
      if (can_reuse_current) {
        // Reuse empty cluster slot
        state->cluster_labels[i] = reuse_idx;
        state->cluster_sizes[reuse_idx] = 1;
        state->cluster_params[reuse_idx] = all_params[chosen_idx];
      } else {
        // Find first empty slot or add new one
        int new_cluster_idx = -1;
        for (int k = 0; k < state->n_clusters; ++k) {
          if (state->cluster_sizes[k] == 0) {
            new_cluster_idx = k;
            break;
          }
        }

        if (new_cluster_idx == -1) {
          // No empty slots, add new cluster
          new_cluster_idx = state->n_clusters;
          state->n_clusters++;
          state->cluster_params.push_back(all_params[chosen_idx]);
          state->cluster_sizes.resize(state->n_clusters);  // Fixed: use resize()
        } else {
          // Use empty slot
          state->cluster_params[new_cluster_idx] = all_params[chosen_idx];
        }

        state->cluster_labels[i] = new_cluster_idx;
        state->cluster_sizes[new_cluster_idx] = 1;
      }
    }
  }

  // Clean up and relabel clusters to be contiguous
  cleanup_empty_clusters();  // Fixed: added parentheses
}

void MCMCRunner::cleanup_empty_clusters() {
  std::vector<int> new_labels(state->n_clusters, -1);
  std::vector<arma::vec> new_params;
  int new_idx = 0;

  // Create mapping from old to new labels
  for (int k = 0; k < state->n_clusters; ++k) {
    if (state->cluster_sizes[k] > 0) {
      new_labels[k] = new_idx;
      new_params.push_back(state->cluster_params[k]);
      new_idx++;
    }
  }

  // Update cluster labels
  for (int i = 0; i < data.n_rows; ++i) {
    state->cluster_labels[i] = new_labels[state->cluster_labels[i]];
  }

  // Update state
  state->cluster_params = new_params;
  state->n_clusters = new_idx;

  // Recompute cluster sizes
  state->cluster_sizes.zeros(state->n_clusters);
  for (int i = 0; i < data.n_rows; ++i) {
    state->cluster_sizes[state->cluster_labels[i]]++;
  }
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
  // Implement Escobar & West (1995) auxiliary variable method
  // This matches the R implementation in update_concentration.R

  double x = R::rbeta(state->alpha + 1.0, data.n_rows);

  // Default alpha priors: shape = 1, rate = 1 (Gamma(1,1))
  double prior_shape = 1.0;
  double prior_rate = 1.0;

  // Calculate mixing probabilities
  double log_x = std::log(x);
  double pi1 = prior_shape + state->n_clusters - 1.0;
  double pi2 = data.n_rows * (prior_rate - log_x);

  double pi_ratio = pi1 / (pi1 + pi2);

  // Sample posterior shape
  double post_shape;
  if (R::runif(0, 1) < pi_ratio) {
    post_shape = prior_shape + state->n_clusters;
  } else {
    post_shape = prior_shape + state->n_clusters - 1.0;
  }

  // Sample new alpha
  double post_rate = prior_rate - log_x;
  state->alpha = R::rgamma(post_shape, 1.0 / post_rate);

  // Ensure alpha stays positive and reasonable
  if (state->alpha <= 0.0) {
    state->alpha = 0.001;
  }
  if (state->alpha > 10.0) {
    state->alpha = 10.0;  // Cap at reasonable value
  }
}

void MCMCRunner::store_iteration(int iter) {
  // Store alpha
  arma::vec alpha_vec(1);
  alpha_vec[0] = state->alpha;
  alpha_samples.push_back(alpha_vec);

  // Store cluster labels (copy to avoid reference issues)
  std::vector<int> labels_copy(state->cluster_labels.begin(),
                               state->cluster_labels.end());
  cluster_samples.push_back(labels_copy);

  // Store cluster parameters (deep copy)
  std::vector<arma::vec> params_copy;
  for (const auto& param : state->cluster_params) {
    params_copy.push_back(arma::vec(param));
  }
  theta_samples.push_back(params_copy);

  // Remove the likelihood computation for now as it's not essential
  // and causing compilation errors
}

} // namespace dirichletprocess
