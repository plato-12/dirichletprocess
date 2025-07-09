# C++ Implementation Function List (Detailed)

This document outlines the key C++ classes and functions planned for the reimplementation of `dirichletprocess` core components. Signatures use `arma` for matrices/vectors and `Rcpp::List` for complex R object interchange.

## 1. Core Data Structures

### 1.1. `MixingDistribution` (Abstract Base Class / Concept)
Represents the component distribution $F(\cdot | \theta)$ and its prior $G_0$.
```cpp
class MixingDistribution {
public:
  virtual ~MixingDistribution() {}

  // Calculate likelihood of data point(s) x given parameters theta
  // x can be arma::vec for univariate, arma::mat for multivariate (rows are obs)
  // theta is Rcpp::List {param1_name: arma::vec/mat, ...}
  virtual arma::vec likelihood(const Rcpp::RObject& x_r, const Rcpp::List& theta) const = 0;

  // Calculate log likelihood (often more stable)
  virtual arma::vec log_likelihood(const Rcpp::RObject& x_r, const Rcpp::List& theta) const = 0;

  // Draw n samples from the prior G0
  // Returns Rcpp::List of parameter arrays {param1_name: arma::vec/mat, ...}
  virtual Rcpp::List prior_draw(int n) const = 0;

  // Draw n samples from the posterior p(theta | data_cluster, G0)
  // data_cluster is arma::vec or arma::mat
  virtual Rcpp::List posterior_draw(const Rcpp::RObject& data_cluster_r, int n_samples) const = 0;

  // Calculate predictive likelihood for new data point(s) x_new (for conjugate models)
  // Returns arma::vec of predictive likelihoods
  virtual arma::vec predictive_likelihood(const Rcpp::RObject& x_new_r) const = 0;

  // For non-conjugate models (used within MetropolisHastings)
  // Calculate prior density p0(theta)
  virtual double prior_density(const Rcpp::List& theta) const = 0; // log prior density often preferred

  // Propose new parameters for Metropolis-Hastings
  virtual Rcpp::List mh_parameter_proposal(const Rcpp::List& current_theta) const = 0;

  // Store prior hyperparameters (e.g., as Rcpp::List or specific members)
  Rcpp::List prior_parameters;
  // Other model-specific fixed parameters (e.g., maxY for Beta, dimension for MVN)
  // double maxY; // Example for Beta
  // int dimension; // Example for MVN
};
```

### 1.2. `DirichletProcess`
Main object holding the DP state.
```cpp
class DirichletProcess {
public:
  // Data
  arma::mat y_data;         // N x D data matrix (D=1 for univariate)
  int n_obs;              // Number of observations
  int data_dim;           // Dimension of data

  // DP Parameters
  double alpha;           // Concentration parameter
  Rcpp::List alpha_prior_params; // Parameters for prior on alpha (e.g., shape, rate)

  // Cluster Assignments & Counts
  arma::uvec cluster_labels;   // Length N, assignments for each data point
  arma::uvec points_per_cluster; // Length K, counts for each cluster
  int num_clusters;         // Current number of active clusters K

  // Cluster Parameters
  // This will be a list of lists, where each inner list is a parameter set {name: value}
  // e.g., params_for_cluster[k]["mu"], params_for_cluster[k]["sigma_sq"]
  Rcpp::List cluster_parameters; // List of Rcpp::List, length K

  // Mixing Distribution (polymorphic)
  std::unique_ptr<MixingDistribution> mixing_dist;

  // For non-conjugate models
  int m_aux_params;         // Number of auxiliary parameters for Neal's Alg 8
  int mh_draws;             // Number of MH draws for parameter updates

  // For conjugate models (precomputed)
  arma::vec predictive_array; // Stores predictive likelihood for each data point

  // Constructor, methods, etc.
  DirichletProcess(const Rcpp::List& r_dp_obj); // Constructor from R list
  Rcpp::List to_R_list() const; // Convert back to R list

  // Core MCMC update steps (will call C++ functions below)
  void update_cluster_components();
  void update_cluster_parameters();
  void update_concentration_parameter();
  void precompute_predictive_array_if_conjugate(); // For conjugate models
};
```

## 2. Likelihood Functions (as methods of specific MixingDistribution subclasses or standalone)

These would ideally be methods of concrete `MixingDistribution` subclasses (e.g., `GaussianMixingDistribution : public MixingDistribution`).
Signatures shown here are for conceptual C++ functions. Input `x_r` is `Rcpp::RObject` to handle `arma::vec` (univariate) or `arma::mat` (multivariate, rows as obs).
`theta` is `Rcpp::List` like `Rcpp::List::create(Rcpp::Named("mu")=mu_val, Rcpp::Named("sigma_sq")=sigma_sq_val)`.

```cpp
// Gaussian: F(y | mu, sigma_sq)
arma::vec likelihood_gaussian(const arma::vec& y, double mu, double sigma_sq);
arma::vec log_likelihood_gaussian(const arma::vec& y, double mu, double sigma_sq);

// Beta: F(y | mu_y, tau, maxY) (using mean/precision on [0, maxY])
arma::vec likelihood_beta(const arma::vec& y, double mu_y, double tau, double maxY);
arma::vec log_likelihood_beta(const arma::vec& y, double mu_y, double tau, double maxY);

// Weibull: F(y | k_shape, lambda_scale)
arma::vec likelihood_weibull(const arma::vec& y, double k_shape, double lambda_scale);
arma::vec log_likelihood_weibull(const arma::vec& y, double k_shape, double lambda_scale);

// Exponential: F(y | lambda_rate)
arma::vec likelihood_exponential(const arma::vec& y, double lambda_rate);
arma::vec log_likelihood_exponential(const arma::vec& y, double lambda_rate);

// Multivariate Normal: F(y_vec | mu_vec, Sigma_mat)
arma::vec likelihood_mvnormal(const arma::mat& y_mat, const arma::vec& mu_vec, const arma::mat& Sigma_mat);
arma::vec log_likelihood_mvnormal(const arma::mat& y_mat, const arma::vec& mu_vec, const arma::mat& Sigma_mat);
```

## 3. Predictive Likelihood Functions (for Conjugate Models)

Calculate $\int F(y_{new} | \phi) dG_0(\phi)$. `prior_params` from `MixingDistribution` object.

```cpp
// Gaussian (Normal-InvGamma G0)
double predictive_gaussian(double y_new, const Rcpp::List& prior_params);

// Exponential (Gamma G0)
double predictive_exponential(double y_new, const Rcpp::List& prior_params);

// Multivariate Normal (Normal-Wishart G0)
double predictive_mvnormal(const arma::vec& y_new_vec, const Rcpp::List& prior_params, int dim);
```

## 4. Prior and Posterior Draw Functions

Return `Rcpp::List` of named parameters, e.g., `Rcpp::List::create(Rcpp::Named("mu")=mu_draws, Rcpp::Named("sigma_sq")=sigma_sq_draws)` where `mu_draws` is `arma::vec` of length `n_samples`.

### 4.1. Prior Draw Functions ($G_0$)
`prior_params` from `MixingDistribution` object.
```cpp
// Gaussian (Normal-InvGamma G0)
Rcpp::List prior_draw_gaussian(const Rcpp::List& prior_params, int n_samples);

// Beta (e.g., U for mu_y, Gamma for tau)
Rcpp::List prior_draw_beta(const Rcpp::List& prior_params, double maxY, int n_samples);

// Weibull (e.g., U for k, InvGamma for lambda)
Rcpp::List prior_draw_weibull(const Rcpp::List& prior_params, int n_samples);

// Exponential (Gamma G0 for rate)
Rcpp::List prior_draw_exponential(const Rcpp::List& prior_params, int n_samples);

// Multivariate Normal (Normal-Wishart G0)
Rcpp::List prior_draw_mvnormal(const Rcpp::List& prior_params, int dim, int n_samples);
```

### 4.2. Posterior Draw Functions (Conjugate Models)
`prior_params` from `MixingDistribution`, `data_cluster` is `arma::vec` or `arma::mat`.
```cpp
// Gaussian (Normal-InvGamma G0)
Rcpp::List posterior_draw_gaussian(const arma::vec& data_cluster, const Rcpp::List& prior_params, int n_samples);

// Exponential (Gamma G0)
Rcpp::List posterior_draw_exponential(const arma::vec& data_cluster, const Rcpp::List& prior_params, int n_samples);

// Multivariate Normal (Normal-Wishart G0)
Rcpp::List posterior_draw_mvnormal(const arma::mat& data_cluster, const Rcpp::List& prior_params, int dim, int n_samples);
```

## 5. MCMC Core Algorithm Components

These functions will likely take a C++ `DirichletProcess& dp_obj_cpp` by reference to modify it, or be methods of the `DirichletProcess` class.
For R interface, they might be wrapped to accept/return `Rcpp::List`.

```cpp
// Updates cluster assignments (Neal's Algorithm 3 for conjugate, 8 for non-conjugate)
void cluster_component_update(DirichletProcess& dp_obj_cpp);
// R-callable wrapper might look like:
// Rcpp::List rcpp_cluster_component_update(Rcpp::List r_dp_obj);

// Updates parameters for each cluster
void cluster_parameter_update(DirichletProcess& dp_obj_cpp);
// R-callable wrapper:
// Rcpp::List rcpp_cluster_parameter_update(Rcpp::List r_dp_obj);

// Updates concentration parameter alpha
void update_alpha_parameter(DirichletProcess& dp_obj_cpp);
// R-callable wrapper (or direct if simpler):
// double rcpp_update_alpha(double current_alpha, int n_obs, int num_clusters, const arma::vec& alpha_prior_params);

// Utility for ClusterComponentUpdate: manages changes to cluster labels, counts, params
void cluster_label_change(DirichletProcess& dp_obj_cpp, 
                          int data_point_idx, 
                          int new_cluster_assignment_logical_idx, // 1 to K (existing) or K+1 to K+M (new/aux)
                          int old_cluster_label, 
                          const Rcpp::List& params_for_new_cluster_if_any); // Parameters if new from aux or posterior draw

// For non-conjugate ClusterParameterUpdate:
Rcpp::List metropolis_hastings_sampler(
    const MixingDistribution& mix_dist,         // To access likelihood, prior_density, mh_proposal
    const Rcpp::RObject& data_cluster_r,        // Data for this cluster
    const Rcpp::List& initial_params,           // Starting parameters for MH
    int num_draws,                              // Number of MH iterations
    const Rcpp::List& mh_proposal_control_params // e.g., step sizes
);

// Helper for MH: Propose new parameters (example for Beta, others needed)
Rcpp::List mh_proposal_beta(const Rcpp::List& current_params, const Rcpp::List& control_params, double maxY);

// Helper for MH: Calculate log prior density (example for Beta, others needed)
double log_prior_density_beta(const Rcpp::List& params, const Rcpp::List& prior_hyperparams, double maxY);
```

## 6. Hierarchical Model Components (Future Scope)

```cpp
// void global_parameter_update_hdp(HDPObject& hdp_obj);
// void update_g0_hdp(HDPObject& hdp_obj);
// void update_gamma_hdp(HDPObject& hdp_obj);
```

## 7. Rcpp Exported Functions (Example Wrappers)

These functions will be exposed to R and will typically convert R objects (Lists) to C++ objects, call the core C++ logic, and convert results back.

```cpp
// [[Rcpp::export]]
Rcpp::List Fit_cpp(Rcpp::List r_dp_obj, int num_iterations, bool progress_bar) {
  // 1. Convert r_dp_obj to C++ DirichletProcess object
  // 2. Loop num_iterations:
  //      call C++ update_cluster_components
  //      call C++ update_cluster_parameters
  //      call C++ update_alpha_parameter
  // 3. Convert C++ DirichletProcess object back to Rcpp::List
  // return updated_r_dp_obj;
  return R_NilValue; // Placeholder
}

// Individual component wrappers if needed for direct R access/testing
// [[Rcpp::export]]
Rcpp::List ClusterComponentUpdate_cpp(Rcpp::List r_dp_obj) { /* ... */ return R_NilValue; }
// [[Rcpp::export]]
Rcpp::List ClusterParameterUpdate_cpp(Rcpp::List r_dp_obj) { /* ... */ return R_NilValue; }
// [[Rcpp::export]]
Rcpp::List UpdateAlpha_cpp(Rcpp::List r_dp_obj) { /* ... */ return R_NilValue; }
```

