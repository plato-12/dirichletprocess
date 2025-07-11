// src/mcmc_runner_manual.cpp
#include "../inst/include/mcmc_runner_manual.h"
#include "../inst/include/utilities.h"
#include <cmath>
#include <algorithm>

namespace dirichletprocess {

void MCMCRunnerManual::initialize_manual_storage() {
  // Pre-allocate storage
  int n_stored = (n_iter - n_burn) / thin;
  alpha_samples.reserve(n_stored);
  cluster_samples.reserve(n_stored);
  theta_samples.reserve(n_stored);
  likelihood_samples.reserve(n_stored);

  // Additional diagnostic storage
  log_posterior_chain.reserve(n_iter);
  entropy_chain.reserve(n_iter);
}

void MCMCRunnerManual::step_cluster_assignments() {
  if (!update_clusters_flag) return;

  // Update auxiliary parameters for Algorithm 8
  update_auxiliary_parameters();

  // Call parent class method with temperature adjustment
  if (temperature != 1.0) {
    // Implement tempered sampling
    for (size_t i = 0; i < data.n_rows; ++i) {
      arma::vec obs = data.row(i).t();
      int current_cluster = state->cluster_labels[i];

      // Remove from current cluster
      if (current_cluster < static_cast<int>(state->cluster_sizes.n_elem)) {
        state->cluster_sizes[current_cluster]--;
      }

      // Calculate tempered probabilities
      std::vector<double> probs;

      // Existing clusters
      for (int k = 0; k < state->n_clusters; ++k) {
        double log_lik = mixing_dist->log_likelihood(obs, state->cluster_params[k]);
        double weight = (k == current_cluster && state->cluster_sizes[k] == 0) ?
        state->alpha / m_auxiliary : state->cluster_sizes[k];
        probs.push_back(weight * std::exp(log_lik / temperature));
      }

      // Auxiliary components
      for (int j = 0; j < m_auxiliary; ++j) {
        double log_lik = mixing_dist->log_likelihood(obs, auxiliary_params[j]);
        probs.push_back((state->alpha / m_auxiliary) * std::exp(log_lik / temperature));
      }

      // Sample new assignment
      int new_assignment = sample_categorical(probs);

      if (new_assignment < state->n_clusters) {
        state->cluster_labels[i] = new_assignment;
        state->cluster_sizes[new_assignment]++;
      } else {
        // New cluster from auxiliary
        int aux_idx = new_assignment - state->n_clusters;
        int new_cluster_idx = state->n_clusters;

        state->cluster_params.push_back(auxiliary_params[aux_idx]);
        state->cluster_labels[i] = new_cluster_idx;

        arma::vec new_sizes(state->cluster_sizes.n_elem + 1);
        new_sizes.head(state->cluster_sizes.n_elem) = state->cluster_sizes;
        new_sizes(new_cluster_idx) = 1;
        state->cluster_sizes = new_sizes;

        state->n_clusters++;
      }
    }

    cleanup_empty_clusters();
  } else {
    // Standard update
    update_cluster_assignments_algorithm8();
  }

  current_iteration++;
}

void MCMCRunnerManual::step_cluster_parameters() {
  if (!update_params_flag) return;

  for (int k = 0; k < state->n_clusters; ++k) {
    if (state->cluster_sizes[k] > 0) {
      // Get cluster data
      std::vector<int> cluster_indices;
      for (size_t i = 0; i < data.n_rows; ++i) {
        if (state->cluster_labels[i] == k) {
          cluster_indices.push_back(i);
        }
      }

      arma::mat cluster_data(cluster_indices.size(), data.n_cols);
      for (size_t idx = 0; idx < cluster_indices.size(); ++idx) {
        cluster_data.row(idx) = data.row(cluster_indices[idx]);
      }

      // Draw from posterior
      arma::vec new_params = mixing_dist->posterior_draw(
        cluster_data, state->cluster_params[k]);

      // Apply bounds if set
      if (param_lower_bounds.n_elem > 0 && check_parameter_bounds(new_params)) {
        state->cluster_params[k] = new_params;
      }
    }
  }
}

void MCMCRunnerManual::step_concentration() {
  if (update_concentration_flag) {
    update_concentration();
  }
}

void MCMCRunnerManual::perform_iteration() {
  step_cluster_assignments();
  step_cluster_parameters();
  step_concentration();

  // Store diagnostics
  log_posterior_chain.push_back(get_log_posterior());
  entropy_chain.push_back(get_clustering_entropy());

  if (current_iteration >= n_burn &&
      (current_iteration - n_burn) % thin == 0) {
    store_iteration(current_iteration);
  }
}

Rcpp::List MCMCRunnerManual::get_current_state() const {
  Rcpp::List params_list;
  for (const auto& param : state->cluster_params) {
    params_list.push_back(param);
  }

  double log_lik = 0.0;
  for (int i = 0; i < data.n_rows; i++) {
    arma::vec obs = data.row(i).t();
    int label = state->cluster_labels[i];
    if (label < state->cluster_params.size()) {
      log_lik += mixing_dist->log_likelihood(obs, state->cluster_params[label]);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("cluster_labels") = state->cluster_labels,
    Rcpp::Named("cluster_params") = params_list,
    Rcpp::Named("cluster_sizes") = state->cluster_sizes,
    Rcpp::Named("alpha") = state->alpha,
    Rcpp::Named("n_clusters") = state->n_clusters,
    Rcpp::Named("iteration") = current_iteration,
    Rcpp::Named("log_likelihood") = log_lik,
    Rcpp::Named("log_posterior") = get_log_posterior(),
    Rcpp::Named("temperature") = temperature
  );
}

void MCMCRunnerManual::set_cluster_labels(const std::vector<int>& new_labels) {
  if (new_labels.size() != state->cluster_labels.size()) {
    Rcpp::stop("New labels must have same length as data");
  }
  state->cluster_labels = new_labels;
  state->update_cluster_counts();
  cleanup_empty_clusters();
}

void MCMCRunnerManual::set_cluster_params(const Rcpp::List& new_params) {
  state->cluster_params.clear();
  for (int i = 0; i < new_params.size(); i++) {
    arma::vec param = Rcpp::as<arma::vec>(new_params[i]);
    if (check_parameter_bounds(param)) {
      state->cluster_params.push_back(param);
    } else {
      Rcpp::warning("Parameter " + std::to_string(i) + " violates bounds, keeping original");
    }
  }
  state->n_clusters = state->cluster_params.size();
}

// Additional Features Implementation

void MCMCRunnerManual::set_parameter_bounds(const arma::vec& lower, const arma::vec& upper) {
  if (lower.n_elem != upper.n_elem) {
    Rcpp::stop("Lower and upper bounds must have same dimension");
  }
  param_lower_bounds = lower;
  param_upper_bounds = upper;
}

Rcpp::List MCMCRunnerManual::get_auxiliary_params() const {
  Rcpp::List aux_list;
  for (const auto& aux : auxiliary_params) {
    aux_list.push_back(aux);
  }
  return aux_list;
}

void MCMCRunnerManual::set_update_flags(bool update_clusters, bool update_params, bool update_alpha) {
  update_clusters_flag = update_clusters;
  update_params_flag = update_params;
  update_concentration_flag = update_alpha;
}

arma::vec MCMCRunnerManual::get_cluster_likelihoods() const {
  arma::vec likes(state->n_clusters);
  for (int k = 0; k < state->n_clusters; k++) {
    double sum_lik = 0.0;
    int count = 0;
    for (int i = 0; i < data.n_rows; i++) {
      if (state->cluster_labels[i] == k) {
        sum_lik += mixing_dist->log_likelihood(
          data.row(i).t(),
          state->cluster_params[k]
        );
        count++;
      }
    }
    likes[k] = count > 0 ? sum_lik : -std::numeric_limits<double>::infinity();
  }
  return likes;
}

arma::mat MCMCRunnerManual::get_cluster_membership_matrix() const {
  arma::mat membership(data.n_rows, state->n_clusters, arma::fill::zeros);
  for (int i = 0; i < data.n_rows; i++) {
    if (state->cluster_labels[i] < state->n_clusters) {
      membership(i, state->cluster_labels[i]) = 1.0;
    }
  }
  return membership;
}

Rcpp::List MCMCRunnerManual::get_cluster_statistics() const {
  Rcpp::List stats;

  for (int k = 0; k < state->n_clusters; k++) {
    // Collect data points in cluster
    std::vector<int> indices;
    for (int i = 0; i < data.n_rows; i++) {
      if (state->cluster_labels[i] == k) {
        indices.push_back(i);
      }
    }

    if (indices.size() > 0) {
      arma::mat cluster_data(indices.size(), data.n_cols);
      for (size_t idx = 0; idx < indices.size(); ++idx) {
        cluster_data.row(idx) = data.row(indices[idx]);
      }

      // Calculate statistics
      arma::vec mean = arma::mean(cluster_data, 0).t();
      arma::mat cov = arma::cov(cluster_data);

      stats.push_back(Rcpp::List::create(
          Rcpp::Named("cluster_id") = k + 1,  // R uses 1-based indexing
          Rcpp::Named("size") = indices.size(),
          Rcpp::Named("mean") = mean,
          Rcpp::Named("covariance") = cov,
          Rcpp::Named("parameters") = state->cluster_params[k],
                                                           Rcpp::Named("likelihood") = get_cluster_likelihoods()[k]
      ));
    }
  }

  return stats;
}

void MCMCRunnerManual::merge_clusters(int cluster1, int cluster2) {
  if (cluster1 >= state->n_clusters || cluster2 >= state->n_clusters) {
    Rcpp::stop("Invalid cluster indices for merging");
  }

  if (cluster1 == cluster2) return;

  // Merge cluster2 into cluster1
  for (int i = 0; i < data.n_rows; i++) {
    if (state->cluster_labels[i] == cluster2) {
      state->cluster_labels[i] = cluster1;
    }
  }

  // Update cluster counts
  state->update_cluster_counts();

  // Re-estimate parameters for merged cluster
  std::vector<int> merged_indices;
  for (int i = 0; i < data.n_rows; i++) {
    if (state->cluster_labels[i] == cluster1) {
      merged_indices.push_back(i);
    }
  }

  if (merged_indices.size() > 0) {
    arma::mat merged_data(merged_indices.size(), data.n_cols);
    for (size_t idx = 0; idx < merged_indices.size(); ++idx) {
      merged_data.row(idx) = data.row(merged_indices[idx]);
    }

    state->cluster_params[cluster1] = mixing_dist->posterior_draw(
      merged_data, state->cluster_params[cluster1]);
  }

  // Clean up empty clusters
  cleanup_empty_clusters();
}

void MCMCRunnerManual::split_cluster(int cluster_id, double split_prob) {
  if (cluster_id >= state->n_clusters) {
    Rcpp::stop("Invalid cluster index for splitting");
  }

  if (state->cluster_sizes[cluster_id] < 2) {
    Rcpp::warning("Cannot split cluster with less than 2 observations");
    return;
  }

  // Create new cluster
  int new_cluster_id = state->n_clusters;
  state->cluster_params.push_back(mixing_dist->prior_draw());

  // Randomly split observations
  for (int i = 0; i < data.n_rows; i++) {
    if (state->cluster_labels[i] == cluster_id) {
      if (R::runif(0, 1) < split_prob) {
        state->cluster_labels[i] = new_cluster_id;
      }
    }
  }

  // Update state
  state->n_clusters++;
  state->update_cluster_counts();

  // Re-estimate parameters for both clusters
  step_cluster_parameters();
}

void MCMCRunnerManual::set_temperature(double temp) {
  if (temp <= 0) {
    Rcpp::stop("Temperature must be positive");
  }
  temperature = temp;
}

void MCMCRunnerManual::set_auxiliary_parameter_count(int m) {
  if (m <= 0) {
    Rcpp::stop("Number of auxiliary parameters must be positive");
  }
  m_auxiliary = m;
  auxiliary_params.resize(m);
  update_auxiliary_parameters();
}

Rcpp::List MCMCRunnerManual::sample_posterior_predictive(int n_samples) {
  Rcpp::List samples;

  for (int i = 0; i < n_samples; i++) {
    // Sample cluster with Chinese Restaurant Process
    std::vector<double> probs;
    for (int k = 0; k < state->n_clusters; k++) {
      probs.push_back(state->cluster_sizes[k]);
    }
    probs.push_back(state->alpha);

    int chosen_cluster = sample_categorical(probs);

    arma::vec params;
    if (chosen_cluster < state->n_clusters) {
      params = state->cluster_params[chosen_cluster];
    } else {
      params = mixing_dist->prior_draw();
    }

    // Sample from likelihood - this would need to be added to MixingDistribution
    // For now, return the parameters
    samples.push_back(params);
  }

  return samples;
}

double MCMCRunnerManual::get_log_posterior() const {
  double log_post = 0.0;

  // Likelihood term
  for (int i = 0; i < data.n_rows; i++) {
    arma::vec obs = data.row(i).t();
    int label = state->cluster_labels[i];
    if (label < state->cluster_params.size()) {
      log_post += mixing_dist->log_likelihood(obs, state->cluster_params[label]);
    }
  }

  // Prior on cluster assignments (CRP)
  for (int k = 0; k < state->n_clusters; k++) {
    if (state->cluster_sizes[k] > 0) {
      log_post += std::log(state->cluster_sizes[k]);
    }
  }
  log_post += state->n_clusters * std::log(state->alpha);

  // Prior on alpha
  log_post += (alpha_prior_shape - 1) * std::log(state->alpha) -
    alpha_prior_rate * state->alpha;

  return log_post;
}

arma::vec MCMCRunnerManual::get_cluster_entropies() const {
  arma::vec entropies(state->n_clusters);

  for (int k = 0; k < state->n_clusters; k++) {
    if (state->cluster_sizes[k] > 0) {
      double p = state->cluster_sizes[k] / (double)data.n_rows;
      entropies[k] = -p * std::log(p);
    } else {
      entropies[k] = 0.0;
    }
  }

  return entropies;
}

double MCMCRunnerManual::get_clustering_entropy() const {
  double entropy = 0.0;

  for (int k = 0; k < state->n_clusters; k++) {
    if (state->cluster_sizes[k] > 0) {
      double p = state->cluster_sizes[k] / (double)data.n_rows;
      entropy -= p * std::log(p);
    }
  }

  return entropy;
}

Rcpp::List MCMCRunnerManual::get_convergence_diagnostics() const {
  // Calculate running statistics
  int n = log_posterior_chain.size();

  if (n < 100) {
    return Rcpp::List::create(
      Rcpp::Named("message") = "Not enough iterations for convergence diagnostics"
    );
  }

  // Split chain for Gelman-Rubin diagnostic
  int split_point = n / 2;

  // First half statistics
  double mean1 = 0.0, var1 = 0.0;
  for (int i = 0; i < split_point; i++) {
    mean1 += log_posterior_chain[i];
  }
  mean1 /= split_point;

  for (int i = 0; i < split_point; i++) {
    var1 += std::pow(log_posterior_chain[i] - mean1, 2);
  }
  var1 /= (split_point - 1);

  // Second half statistics
  double mean2 = 0.0, var2 = 0.0;
  for (int i = split_point; i < n; i++) {
    mean2 += log_posterior_chain[i];
  }
  mean2 /= (n - split_point);

  for (int i = split_point; i < n; i++) {
    var2 += std::pow(log_posterior_chain[i] - mean2, 2);
  }
  var2 /= (n - split_point - 1);

  // Approximate R-hat
  double W = (var1 + var2) / 2.0;
  double B = n * std::pow(mean1 - mean2, 2) / 2.0;
  double var_est = W + B / n;
  double R_hat = std::sqrt(var_est / W);

  // Effective sample size (rough approximation)
  double autocorr = 0.0;
  double overall_mean = (mean1 + mean2) / 2.0;
  for (int i = 1; i < n - 1; i++) {
    autocorr += (log_posterior_chain[i] - overall_mean) *
      (log_posterior_chain[i-1] - overall_mean);
  }
  autocorr /= ((n - 1) * var_est);
  double ess = n / (1 + 2 * autocorr);

  return Rcpp::List::create(
    Rcpp::Named("iterations_completed") = current_iteration,
    Rcpp::Named("log_posterior_mean") = (mean1 + mean2) / 2.0,
    Rcpp::Named("log_posterior_sd") = std::sqrt(var_est),
    Rcpp::Named("R_hat") = R_hat,
    Rcpp::Named("effective_sample_size") = ess,
    Rcpp::Named("acceptance_rate") = 1.0, // Would need to track this
    Rcpp::Named("mean_clusters") = arma::mean(n_clusters_chain),
    Rcpp::Named("mean_entropy") = arma::mean(entropy_chain)
  );
}

Rcpp::List MCMCRunnerManual::get_results() const {
  int n_samples = alpha_samples.size();

  Rcpp::NumericMatrix labels_matrix(n_samples, data.n_rows);
  for (int i = 0; i < n_samples; i++) {
    for (int j = 0; j < data.n_rows; j++) {
      labels_matrix(i, j) = cluster_samples[i][j] + 1;
    }
  }

  Rcpp::List theta_list;
  for (int i = 0; i < n_samples; i++) {
    Rcpp::List iter_params;
    for (size_t j = 0; j < theta_samples[i].size(); j++) {
      iter_params.push_back(theta_samples[i][j]);
    }
    theta_list.push_back(iter_params);
  }

  return Rcpp::List::create(
    Rcpp::Named("cluster_labels") = state->cluster_labels,
    Rcpp::Named("cluster_params") = state->cluster_params,
    Rcpp::Named("alpha") = state->alpha,
    Rcpp::Named("n_clusters") = state->n_clusters,
    Rcpp::Named("labels_chain") = labels_matrix,
    Rcpp::Named("alpha_chain") = alpha_samples,
    Rcpp::Named("theta_chain") = theta_list,
    Rcpp::Named("n_clusters_chain") = n_clusters_chain,
    Rcpp::Named("likelihood_chain") = likelihood_samples,
    Rcpp::Named("log_posterior_chain") = log_posterior_chain,
    Rcpp::Named("entropy_chain") = entropy_chain,
    Rcpp::Named("iterations_completed") = current_iteration,
    Rcpp::Named("convergence_diagnostics") = get_convergence_diagnostics()
  );
}

void MCMCRunnerManual::update_auxiliary_parameters() {
  for (int j = 0; j < m_auxiliary; j++) {
    auxiliary_params[j] = mixing_dist->prior_draw();
  }
}

bool MCMCRunnerManual::check_parameter_bounds(const arma::vec& params) const {
  if (param_lower_bounds.n_elem == 0) return true;

  if (params.n_elem != param_lower_bounds.n_elem) {
    Rcpp::warning("Parameter dimension mismatch with bounds");
    return false;
  }

  for (size_t i = 0; i < params.n_elem; i++) {
    if (params[i] < param_lower_bounds[i] || params[i] > param_upper_bounds[i]) {
      return false;
    }
  }

  return true;
}

} // namespace dirichletprocess
