#include <RcppArmadillo.h>
#include "../inst/include/mcmc_runner.h"
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/gaussian_mixing.h"
#include <algorithm>

namespace dirichletprocess {

// Helper function to sample from categorical distribution
int sample_categorical(const arma::vec& probs) {
  arma::vec cumprobs = arma::cumsum(probs / arma::sum(probs));
  double u = R::runif(0, 1);

  for (arma::uword i = 0; i < cumprobs.n_elem; ++i) {
    if (u <= cumprobs[i]) {
      return static_cast<int>(i);
    }
  }
  return static_cast<int>(probs.n_elem - 1);
}

MCMCRunner::MCMCRunner(const arma::mat& data,
                       const Rcpp::List& mixing_dist_params,
                       const Rcpp::List& mcmc_params)
  : data(data) {

  // Validate inputs
  if (data.n_rows == 0) {
    Rcpp::stop("Data cannot be empty");
  }

  if (data.has_nan() || data.has_inf()) {
    Rcpp::stop("Data contains NA or Inf values");
  }

  // Extract MCMC parameters
  n_iter = Rcpp::as<int>(mcmc_params["n_iter"]);
  n_burn = Rcpp::as<int>(mcmc_params["n_burn"]);
  thin = Rcpp::as<int>(mcmc_params["thin"]);
  update_concentration_flag = Rcpp::as<bool>(mcmc_params["update_concentration"]);

  // Validate MCMC parameters
  if (n_iter <= 0) {
    Rcpp::stop("n_iter must be positive");
  }
  if (n_burn >= n_iter) {
    Rcpp::stop("n_burn must be less than n_iter");
  }
  if (thin <= 0) {
    Rcpp::stop("thin must be positive");
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

  // Pre-allocate storage for ALL iterations (not just post-burn)
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
    // Update cluster assignments (Chinese Restaurant Process)
    update_cluster_assignments();

    // Update cluster parameters
    update_cluster_parameters();

    // Update concentration parameter
    if (update_concentration_flag) {
      update_concentration();
    }

    // Store current iteration (store ALL iterations, not just post-burn)
    store_iteration(iter);
  }

  // Convert stored samples to proper format for R
  // Extract only post-burn samples according to thinning
  int n_stored = 0;
  std::vector<arma::vec> alpha_chain;
  std::vector<std::vector<int>> labels_chain;
  std::vector<std::vector<arma::vec>> theta_chain;
  std::vector<int> n_clusters_chain;

  for (int iter = n_burn; iter < n_iter; iter += thin) {
    alpha_chain.push_back(alpha_samples[iter]);
    labels_chain.push_back(cluster_samples[iter]);
    theta_chain.push_back(theta_samples[iter]);

    // Count unique clusters
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
    for (int j = 0; j < data.n_rows; ++j) {
      // Convert to 1-indexed for R
      labels_matrix(i, j) = labels_chain[i][j] + 1;
    }
  }

  // Get final state
  int final_iter = n_iter - 1;
  arma::vec final_labels = arma::conv_to<arma::vec>::from(cluster_samples[final_iter]) + 1;

  // Return results with proper structure
  return Rcpp::List::create(
    Rcpp::Named("cluster_labels") = labels_matrix,      // Matrix of labels (each row is an iteration)
    Rcpp::Named("alpha") = alpha_vector,                // Vector of alpha values
    Rcpp::Named("theta") = theta_chain,                 // List of cluster parameters
    Rcpp::Named("n_clusters") = n_clusters_chain,       // Vector of cluster counts
    Rcpp::Named("final_labels") = final_labels,         // Final cluster assignments
    Rcpp::Named("final_n_clusters") = state->n_clusters // Final number of clusters
  );
}

void MCMCRunner::update_cluster_assignments() {
  int n = data.n_rows;

  for (int i = 0; i < n; ++i) {
    arma::vec obs = data.row(i).t();
    int current_cluster = state->cluster_labels[i];

    // Remove observation from current cluster
    state->cluster_sizes[current_cluster]--;

    // Handle empty cluster
    if (state->cluster_sizes[current_cluster] == 0) {
      // Remove empty cluster
      state->cluster_params.erase(state->cluster_params.begin() + current_cluster);
      state->cluster_sizes.shed_row(current_cluster);

      // Relabel clusters
      for (int j = 0; j < n; ++j) {
        if (state->cluster_labels[j] > current_cluster) {
          state->cluster_labels[j]--;
        }
      }
      state->n_clusters--;
    }

    // Chinese Restaurant Process
    arma::vec probs(state->n_clusters + 1);

    // Probability of joining existing clusters
    for (int k = 0; k < state->n_clusters; ++k) {
      double log_lik = mixing_dist->log_likelihood(obs, state->cluster_params[k]);
      probs[k] = state->cluster_sizes[k] * std::exp(log_lik);
    }

    // Probability of creating new cluster
    // Compute predictive likelihood by integrating over prior
    arma::vec prior_sample = mixing_dist->prior_draw();
    double log_lik_new = mixing_dist->log_likelihood(obs, prior_sample);
    probs[state->n_clusters] = state->alpha * std::exp(log_lik_new);

    // Sample new cluster assignment
    int new_cluster = sample_categorical(probs);

    if (new_cluster == state->n_clusters) {
      // Create new cluster
      state->cluster_params.push_back(prior_sample);
      state->cluster_sizes.resize(state->n_clusters + 1);
      state->cluster_sizes[state->n_clusters] = 1;
      state->cluster_labels[i] = state->n_clusters;
      state->n_clusters++;
    } else {
      // Join existing cluster
      state->cluster_labels[i] = new_cluster;
      state->cluster_sizes[new_cluster]++;
    }
  }

  // Update cluster count
  state->update_cluster_counts();
}

void MCMCRunner::update_cluster_parameters() {
  for (int k = 0; k < state->n_clusters; ++k) {
    // Get data points in cluster k
    arma::mat cluster_data;
    int count = 0;

    for (int i = 0; i < data.n_rows; ++i) {
      if (state->cluster_labels[i] == k) {
        count++;
      }
    }

    if (count > 0) {
      cluster_data.set_size(count, data.n_cols);
      int idx = 0;
      for (int i = 0; i < data.n_rows; ++i) {
        if (state->cluster_labels[i] == k) {
          cluster_data.row(idx++) = data.row(i);
        }
      }

      // Draw from posterior
      arma::vec prior_params = mixing_dist->prior_draw();
      state->cluster_params[k] = mixing_dist->posterior_draw(cluster_data, prior_params);
    }
  }
}

void MCMCRunner::update_concentration() {
  // Auxiliary variable method for updating alpha
  double a = 1.0;  // Gamma prior shape
  double b = 1.0;  // Gamma prior rate

  // Auxiliary variable method (Escobar & West 1995)
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
  // Store alpha
  alpha_samples.push_back(arma::vec{state->alpha});

  // Store cluster labels
  std::vector<int> labels(state->cluster_labels.begin(), state->cluster_labels.end());
  cluster_samples.push_back(labels);

  // Store cluster parameters (deep copy)
  std::vector<arma::vec> params_copy;
  for (const auto& param : state->cluster_params) {
    params_copy.push_back(param);
  }
  theta_samples.push_back(params_copy);
}

// DPState implementation
DPState::DPState(int n_obs, double initial_alpha)
  : cluster_labels(n_obs), alpha(initial_alpha), n_clusters(1) {
  cluster_sizes.set_size(1);
  cluster_sizes[0] = n_obs;
  std::fill(cluster_labels.begin(), cluster_labels.end(), 0);
}

void DPState::update_cluster_counts() {
  // Count unique clusters and update sizes
  std::set<int> unique_labels;
  for (arma::uword i = 0; i < cluster_labels.n_elem; ++i) {
    unique_labels.insert(cluster_labels[i]);
  }
  n_clusters = unique_labels.size();

  // Recompute cluster sizes
  cluster_sizes.zeros(n_clusters);
  for (arma::uword i = 0; i < cluster_labels.n_elem; ++i) {
    cluster_sizes[cluster_labels[i]]++;
  }
}

} // namespace dirichletprocess
