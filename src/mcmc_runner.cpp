#include <RcppArmadillo.h>
#include "../inst/include/mcmc_runner.h"
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/gaussian_mixing.h"
#include <algorithm>

namespace dirichletprocess {

int sample_categorical(const arma::vec& probs) {
  // Check for all zeros or negative probabilities
  if (arma::all(probs <= 0)) {
    Rcpp::stop("All probabilities are zero or negative");
  }

  // Normalize probabilities to avoid numerical issues
  arma::vec normalized_probs = probs / arma::sum(probs);

  // Check for NaN after normalization
  if (normalized_probs.has_nan()) {
    Rcpp::stop("Probability normalization resulted in NaN");
  }

  double u = R::runif(0, 1);
  double cumsum = 0.0;

  for (arma::uword i = 0; i < normalized_probs.n_elem; ++i) {
    cumsum += normalized_probs[i];
    if (u <= cumsum) {
      return static_cast<int>(i);
    }
  }

  // This should never happen, but return last index as fallback
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
    for (arma::uword j = 0; j < data.n_rows; j++) {
      // Convert to 1-indexed for R
      labels_matrix(i, j) = labels_chain[i][j] + 1;
    }
  }

  // Get final state
  int final_iter = n_iter - 1;
  arma::vec final_labels = arma::conv_to<arma::vec>::from(cluster_samples[final_iter]) + 1;

  // Extract final cluster parameters
  std::vector<arma::vec> final_theta = theta_samples[final_iter];

  // Convert theta_chain to List for R
  Rcpp::List theta_chain_list(n_stored);
  for (int i = 0; i < n_stored; ++i) {
    Rcpp::List iter_params(theta_chain[i].size());
    for (size_t j = 0; j < theta_chain[i].size(); ++j) {
      iter_params[j] = theta_chain[i][j];
    }
    theta_chain_list[i] = iter_params;
  }

  // Convert final theta to List
  Rcpp::List final_theta_list(final_theta.size());
  for (size_t i = 0; i < final_theta.size(); ++i) {
    final_theta_list[i] = final_theta[i];
  }

  // Return results with proper structure
  return Rcpp::List::create(
    // Primary results
    Rcpp::Named("cluster_labels") = labels_matrix,      // Matrix of labels (each row is an iteration)
    Rcpp::Named("alpha") = alpha_vector,                // Vector of alpha values
    Rcpp::Named("theta") = final_theta_list,            // FINAL cluster parameters (not chain!)
    Rcpp::Named("n_clusters") = n_clusters_chain,       // Vector of cluster counts

    // Additional fields for compatibility
    Rcpp::Named("final_labels") = final_labels,         // Final cluster assignments
    Rcpp::Named("final_n_clusters") = state->n_clusters,// Final number of clusters

    // Chain fields that R code expects
    Rcpp::Named("alpha_chain") = alpha_vector,          // Same as alpha for compatibility
    Rcpp::Named("labels_chain") = labels_matrix,        // Same as cluster_labels
    Rcpp::Named("theta_chain") = theta_chain_list       // Chain of cluster parameters
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

      // Check for invalid log likelihood
      if (!std::isfinite(log_lik)) {
        log_lik = -1e10; // Set to very small value
      }

      probs[k] = state->cluster_sizes[k] * std::exp(log_lik);
    }

    // Probability of creating new cluster
    // Draw multiple auxiliary parameters and average their likelihoods
    const int n_aux = 5; // Number of auxiliary parameters
    double avg_log_lik = 0.0;

    for (int aux = 0; aux < n_aux; ++aux) {
      arma::vec aux_params = mixing_dist->prior_draw();
      double log_lik = mixing_dist->log_likelihood(obs, aux_params);

      if (std::isfinite(log_lik)) {
        avg_log_lik += log_lik;
      } else {
        avg_log_lik += -1e10;
      }
    }
    avg_log_lik /= n_aux;

    probs[state->n_clusters] = state->alpha * std::exp(avg_log_lik);

    // Add small epsilon to avoid all-zero probabilities
    probs += 1e-10;

    // Sample new cluster assignment
    int new_cluster = sample_categorical(probs);

    if (new_cluster == state->n_clusters) {
      // Create new cluster with proper parameter initialization
      arma::vec new_params = mixing_dist->prior_draw();
      state->cluster_params.push_back(new_params);
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

      // Draw from posterior (don't use prior_params here)
      state->cluster_params[k] = mixing_dist->posterior_draw(cluster_data, state->cluster_params[k]);
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

} // namespace dirichletprocess
