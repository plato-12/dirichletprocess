// src/mcmc_runner.cpp
#include "../inst/include/mcmc_runner.h"
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/utilities.h"
#include <set>
#include <algorithm>
#include <numeric>

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

  // Extract alpha prior parameters
  if (mcmc_params.containsElementNamed("alpha_prior_shape")) {
    alpha_prior_shape = Rcpp::as<double>(mcmc_params["alpha_prior_shape"]);
  } else {
    alpha_prior_shape = 1.0;  // Default Gamma(1,1)
  }

  if (mcmc_params.containsElementNamed("alpha_prior_rate")) {
    alpha_prior_rate = Rcpp::as<double>(mcmc_params["alpha_prior_rate"]);
  } else {
    alpha_prior_rate = 1.0;  // Default Gamma(1,1)
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
  likelihood_samples.reserve(n_iter);
}

int MCMCRunner::sample_categorical(const std::vector<double>& probs) {
  double u = R::runif(0, 1);
  double cumsum = 0.0;

  for (size_t i = 0; i < probs.size(); ++i) {
    cumsum += probs[i];
    if (u <= cumsum) {
      return i;
    }
  }

  return probs.size() - 1;  // Fallback to last category
}

void MCMCRunner::cleanup_empty_clusters() {
  std::vector<int> new_labels(state->cluster_labels.size());
  std::vector<arma::vec> new_params;

  int new_idx = 0;
  std::map<int, int> old_to_new;

  // Build mapping and new parameters
  for (int k = 0; k < state->n_clusters; ++k) {
    if (k < static_cast<int>(state->cluster_sizes.n_elem) && state->cluster_sizes[k] > 0) {
      old_to_new[k] = new_idx;
      if (k < static_cast<int>(state->cluster_params.size())) {
        new_params.push_back(state->cluster_params[k]);
      }
      new_idx++;
    }
  }

  // Update labels
  for (size_t i = 0; i < state->cluster_labels.size(); ++i) {
    if (old_to_new.find(state->cluster_labels[i]) != old_to_new.end()) {
      new_labels[i] = old_to_new[state->cluster_labels[i]];
    } else {
      // This shouldn't happen, but handle gracefully
      new_labels[i] = 0;
    }
  }

  // Update state
  state->cluster_labels = new_labels;
  state->cluster_params = new_params;
  state->n_clusters = new_params.size();
  state->update_cluster_counts();
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
  std::vector<double> likelihood_chain;

  for (int iter = n_burn; iter < n_iter; iter += thin) {
    alpha_chain.push_back(alpha_samples[iter]);
    labels_chain.push_back(cluster_samples[iter]);
    theta_chain.push_back(theta_samples[iter]);
    likelihood_chain.push_back(likelihood_samples[iter]);

    std::set<int> unique_labels(cluster_samples[iter].begin(),
                                cluster_samples[iter].end());
    n_clusters_chain.push_back(unique_labels.size());
    n_stored++;
  }

  // Convert to matrices for R
  arma::mat labels_matrix(n_stored, data.n_rows);
  arma::vec alpha_vector(n_stored);
  arma::vec likelihood_vector(n_stored);

  for (int i = 0; i < n_stored; ++i) {
    alpha_vector[i] = alpha_chain[i][0];
    likelihood_vector[i] = likelihood_chain[i];
    for (size_t j = 0; j < data.n_rows; ++j) {
      labels_matrix(i, j) = labels_chain[i][j] + 1;  // Convert to 1-indexed for R
    }
  }

  // Convert theta to list format
  Rcpp::List theta_list(n_stored);
  for (int i = 0; i < n_stored; ++i) {
    Rcpp::List iter_params(theta_chain[i].size());
    for (size_t j = 0; j < theta_chain[i].size(); ++j) {
      iter_params[j] = theta_chain[i][j];
    }
    theta_list[i] = iter_params;
  }

  // Convert n_clusters to vector
  Rcpp::IntegerVector n_clusters_vector(n_clusters_chain.begin(),
                                        n_clusters_chain.end());

  return Rcpp::List::create(
    Rcpp::Named("labels_chain") = labels_matrix,
    Rcpp::Named("alpha_chain") = alpha_vector,
    Rcpp::Named("theta_chain") = theta_list,
    Rcpp::Named("n_clusters") = n_clusters_vector,
    Rcpp::Named("likelihood_chain") = likelihood_vector,
    Rcpp::Named("cluster_labels") = labels_chain,
    Rcpp::Named("alpha") = alpha_chain,
    Rcpp::Named("theta") = theta_chain
  );
}

void MCMCRunner::update_cluster_assignments_algorithm8() {
  cleanup_empty_clusters();  // First clean up any empty clusters

  for (int i = 0; i < data.n_rows; ++i) {
    arma::vec obs = data.row(i).t();
    int current_cluster = state->cluster_labels[i];

    // Remove observation from current cluster
    if (current_cluster < static_cast<int>(state->cluster_sizes.n_elem)) {
      state->cluster_sizes[current_cluster]--;
    }

    // Prepare probabilities for existing clusters and auxiliary parameters
    std::vector<double> probs;
    std::vector<arma::vec> candidate_params;
    std::vector<int> candidate_indices;  // Track if it's existing cluster or new

    // Add existing clusters
    for (int k = 0; k < state->n_clusters; ++k) {
      if (k < static_cast<int>(state->cluster_params.size())) {
        double log_lik = mixing_dist->log_likelihood(obs, state->cluster_params[k]);
        double weight;

        if (k == current_cluster && state->cluster_sizes[k] == 0) {
          // If this cluster is now empty, treat it like an auxiliary
          weight = state->alpha / m_auxiliary;
        } else {
          weight = state->cluster_sizes[k];
        }

        probs.push_back(weight * std::exp(log_lik));
        candidate_params.push_back(state->cluster_params[k]);
        candidate_indices.push_back(k);
      }
    }

    int n_existing = probs.size();

    // Add m auxiliary parameters
    for (int j = 0; j < m_auxiliary; ++j) {
      arma::vec aux_param = mixing_dist->prior_draw();
      double log_lik = mixing_dist->log_likelihood(obs, aux_param);
      probs.push_back((state->alpha / m_auxiliary) * std::exp(log_lik));
      candidate_params.push_back(aux_param);
      candidate_indices.push_back(-1);  // Mark as new cluster
    }

    // Handle numerical issues
    double max_prob = *std::max_element(probs.begin(), probs.end());
    if (max_prob > 1e100) {
      for (auto& p : probs) {
        p /= max_prob;
      }
    }

    // Normalize and sample
    double prob_sum = std::accumulate(probs.begin(), probs.end(), 0.0);
    if (prob_sum <= 0.0 || !std::isfinite(prob_sum)) {
      // Fallback to uniform
      std::fill(probs.begin(), probs.end(), 1.0 / probs.size());
    } else {
      for (auto& p : probs) {
        p /= prob_sum;
      }
    }

    int chosen_idx = sample_categorical(probs);

    // Assign to cluster
    if (chosen_idx < n_existing && candidate_indices[chosen_idx] >= 0) {
      // Existing cluster
      int cluster_idx = candidate_indices[chosen_idx];
      state->cluster_labels[i] = cluster_idx;
      state->cluster_sizes[cluster_idx]++;
    } else {
      // New cluster from auxiliary parameter
      int new_cluster_idx;

      // Check if we can reuse the empty current cluster
      if (current_cluster < state->n_clusters && state->cluster_sizes[current_cluster] == 0) {
        new_cluster_idx = current_cluster;
        state->cluster_params[new_cluster_idx] = candidate_params[chosen_idx];
      } else {
        // Create entirely new cluster
        new_cluster_idx = state->n_clusters;
        state->cluster_params.push_back(candidate_params[chosen_idx]);
        state->cluster_sizes.conservativeResize(state->n_clusters + 1);
        state->n_clusters++;
      }

      state->cluster_labels[i] = new_cluster_idx;
      state->cluster_sizes[new_cluster_idx] = 1;
    }
  }

  cleanup_empty_clusters();  // Clean up after all assignments
}

void MCMCRunner::update_cluster_parameters() {
  for (int k = 0; k < state->n_clusters; ++k) {
    if (state->cluster_sizes[k] > 0) {
      // Get indices of observations in this cluster
      std::vector<int> cluster_indices;
      for (int i = 0; i < data.n_rows; ++i) {
        if (state->cluster_labels[i] == k) {
          cluster_indices.push_back(i);
        }
      }

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
  // Escobar & West (1995) auxiliary variable method
  double x = R::rbeta(state->alpha + 1.0, data.n_rows);

  // Use the proper priors
  double log_x = std::log(x);
  if (!std::isfinite(log_x)) {
    log_x = -10.0;  // Fallback for numerical stability
  }

  double pi1 = alpha_prior_shape + state->n_clusters - 1.0;
  double pi2 = data.n_rows * (alpha_prior_rate - log_x);

  double pi_ratio = pi1 / (pi1 + pi2);
  if (!std::isfinite(pi_ratio)) {
    pi_ratio = 0.5;  // Fallback
  }

  // Sample posterior shape
  double post_shape;
  if (R::runif(0, 1) < pi_ratio) {
    post_shape = alpha_prior_shape + state->n_clusters;
  } else {
    post_shape = alpha_prior_shape + state->n_clusters - 1.0;
  }

  // Sample new alpha
  double post_rate = alpha_prior_rate - log_x;
  if (post_rate <= 0.0) {
    post_rate = 0.001;  // Ensure positive rate
  }

  state->alpha = R::rgamma(post_shape, 1.0 / post_rate);

  // Ensure alpha stays in reasonable range
  if (state->alpha <= 0.0) {
    state->alpha = 0.001;
  }
  if (state->alpha > 100.0) {
    state->alpha = 100.0;  // Cap at reasonable value
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

  // Calculate and store likelihood
  double log_lik = 0.0;
  for (int i = 0; i < data.n_rows; ++i) {
    arma::vec obs = data.row(i).t();
    int cluster = state->cluster_labels[i];
    if (cluster >= 0 && cluster < static_cast<int>(state->cluster_params.size())) {
      log_lik += mixing_dist->log_likelihood(obs, state->cluster_params[cluster]);
    }
  }
  likelihood_samples.push_back(log_lik);
}

} // namespace dirichletprocess
