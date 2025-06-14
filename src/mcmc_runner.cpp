// src/mcmc_runner.cpp
#include "mcmc_runner.h"
#include "gaussian_mixing.h"
#include <algorithm>

namespace dirichletprocess {

MCMCRunner::MCMCRunner(const arma::mat& data,
                       const Rcpp::List& mixing_dist_params,
                       const Rcpp::List& mcmc_params)
  : data(data) {

  // Extract MCMC parameters
  n_iter = mcmc_params["n_iter"];
  n_burn = mcmc_params["n_burn"];
  thin = mcmc_params["thin"];
  update_concentration = mcmc_params["update_concentration"];

  // Create mixing distribution
  std::string dist_type = mixing_dist_params["type"];
  mixing_dist = MixingDistribution::create(dist_type, mixing_dist_params);

  // Initialize state
  double initial_alpha = mcmc_params["alpha"];
  state = std::make_unique<DPState>(data.n_rows, initial_alpha);

  // Pre-allocate storage
  int n_store = (n_iter - n_burn) / thin;
  alpha_samples.reserve(n_store);
  cluster_samples.reserve(n_store);
  theta_samples.reserve(n_store);
}

Rcpp::List MCMCRunner::run() {
  // Initialize with one cluster containing all data
  std::fill(state->cluster_labels.begin(), state->cluster_labels.end(), 0);
  state->update_cluster_counts();

  // Initialize cluster parameters
  state->cluster_params.resize(1);
  state->cluster_params[0] = mixing_dist->prior_draw();

  // Main MCMC loop
  for (int iter = 0; iter < n_iter; ++iter) {
    update_cluster_assignments();
    update_cluster_parameters();

    if (update_concentration) {
      update_concentration();
    }

    // Store samples after burn-in and according to thinning
    if (iter >= n_burn && (iter - n_burn) % thin == 0) {
      store_iteration(iter);
    }
  }

  // Package results
  return Rcpp::List::create(
    Rcpp::Named("cluster_labels") = cluster_samples,
    Rcpp::Named("alpha") = alpha_samples,
    Rcpp::Named("theta") = theta_samples,
    Rcpp::Named("n_clusters") = state->n_clusters
  );
}

void MCMCRunner::update_cluster_assignments() {
  int n = data.n_rows;

  for (int i = 0; i < n; ++i) {
    arma::vec obs = data.row(i).t();
    int current_cluster = state->cluster_labels[i];

    // Remove observation from current cluster
    state->cluster_sizes[current_cluster]--;
    if (state->cluster_sizes[current_cluster] == 0) {
      // Remove empty cluster
      state->cluster_params.erase(state->cluster_params.begin() + current_cluster);
      state->cluster_sizes.shed_row(current_cluster);

      // Adjust labels
      for (int j = 0; j < n; ++j) {
        if (state->cluster_labels[j] > current_cluster) {
          state->cluster_labels[j]--;
        }
      }
      state->n_clusters--;
    }

    // Calculate probabilities for existing clusters and new cluster
    arma::vec probs(state->n_clusters + 1);

    // Existing clusters
    for (int k = 0; k < state->n_clusters; ++k) {
      double log_prob = std::log(state->cluster_sizes[k]) +
        mixing_dist->log_likelihood(obs, state->cluster_params[k]);
      probs[k] = std::exp(log_prob);
    }

    // New cluster
    arma::vec new_params = mixing_dist->prior_draw();
    double log_prob_new = std::log(state->alpha) +
      mixing_dist->log_likelihood(obs, new_params);
    probs[state->n_clusters] = std::exp(log_prob_new);

    // Normalize and sample
    probs = probs / arma::sum(probs);
    int new_cluster = sample_categorical(probs);

    // Assign to cluster
    if (new_cluster == state->n_clusters) {
      // Create new cluster
      state->cluster_params.push_back(new_params);
      state->cluster_sizes.resize(state->n_clusters + 1);
      state->cluster_sizes[state->n_clusters] = 1;
      state->n_clusters++;
    } else {
      state->cluster_sizes[new_cluster]++;
    }

    state->cluster_labels[i] = new_cluster;
  }
}

void MCMCRunner::update_cluster_parameters() {
  for (int k = 0; k < state->n_clusters; ++k) {
    // Get data in cluster k
    arma::uvec cluster_idx = arma::find(state->cluster_labels == k);
    arma::mat cluster_data = data.rows(cluster_idx);

    // Draw from posterior
    state->cluster_params[k] = mixing_dist->posterior_draw(cluster_data,
                                                           state->cluster_params[k]);
  }
}

void MCMCRunner::update_concentration() {
  // Implement concentration parameter update (e.g., using auxiliary variable method)
  // This is a simplified version - you'd want the full auxiliary variable sampler
  double a = 1.0; // hyperparameter
  double b = 1.0; // hyperparameter

  // Auxiliary variable method
  double eta = R::rbeta(state->alpha + 1, data.n_rows);
  double pi_eta = (a + state->n_clusters - 1) /
    (a + state->n_clusters - 1 + data.n_rows * (b - std::log(eta)));

  if (R::runif(0, 1) < pi_eta) {
    state->alpha = R::rgamma(a + state->n_clusters, 1.0 / (b - std::log(eta)));
  } else {
    state->alpha = R::rgamma(a + state->n_clusters - 1, 1.0 / (b - std::log(eta)));
  }
}

// Helper function
int sample_categorical(const arma::vec& probs) {
  double u = R::runif(0, 1);
  double cumsum = 0;
  for (size_t i = 0; i < probs.n_elem; ++i) {
    cumsum += probs[i];
    if (u <= cumsum) return i;
  }
  return probs.n_elem - 1;
}

} // namespace dirichletprocess
