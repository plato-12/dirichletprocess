// src/MVNormalDistribution.cpp
#include "../inst/include/MVNormalDistribution.h"
#include "../inst/include/RcppConversions.h"
#include <RcppArmadillo.h>

namespace dp {

// Helper function to parse covariance model from string
CovarianceModel parseCovarianceModel(const std::string& model) {
  if (model == "E") return CovarianceModel::E;
  else if (model == "V") return CovarianceModel::V;
  else if (model == "EII") return CovarianceModel::EII;
  else if (model == "VII") return CovarianceModel::VII;
  else if (model == "EEI") return CovarianceModel::EEI;
  else if (model == "VEI") return CovarianceModel::VEI;
  else if (model == "EVI") return CovarianceModel::EVI;
  else if (model == "VVI") return CovarianceModel::VVI;
  else return CovarianceModel::FULL;
}

// MVNormalMixingDistribution implementation
MVNormalMixingDistribution::MVNormalMixingDistribution(const Rcpp::List& priorParams) {
  distribution = "mvnormal";
  conjugate = true;
  priorParameters = priorParams;

  // Extract prior parameters
  if (priorParams.containsElementNamed("mu0")) {
    Rcpp::NumericVector mu0_vec = Rcpp::as<Rcpp::NumericVector>(priorParams["mu0"]);
    mu0 = arma::vec(mu0_vec.begin(), mu0_vec.size());
  }

  if (priorParams.containsElementNamed("kappa0")) {
    kappa0 = Rcpp::as<double>(priorParams["kappa0"]);
  }

  if (priorParams.containsElementNamed("Lambda")) {
    Lambda = Rcpp::as<arma::mat>(priorParams["Lambda"]);
  }

  if (priorParams.containsElementNamed("nu")) {
    nu = Rcpp::as<double>(priorParams["nu"]);
  }

  // Extract covariance model
  if (priorParams.containsElementNamed("covModel")) {
    std::string modelStr = Rcpp::as<std::string>(priorParams["covModel"]);
    covModel = parseCovarianceModel(modelStr);
  } else {
    covModel = CovarianceModel::FULL;
  }
}

MVNormalMixingDistribution::~MVNormalMixingDistribution() {
  // Destructor
}

// Get number of covariance parameters for the model
int MVNormalMixingDistribution::getNumCovParams(int d) const {
  switch (covModel) {
  case CovarianceModel::E:
    return 1;  // One variance parameter
  case CovarianceModel::V:
    return 1;  // One variance parameter per observation
  case CovarianceModel::EII:
    return 1;  // One volume parameter
  case CovarianceModel::VII:
    return 1;  // One volume parameter per cluster
  case CovarianceModel::EEI:
    return d;  // Diagonal elements (same across clusters)
  case CovarianceModel::VEI:
    return d + 1;  // Volume + diagonal shape
  case CovarianceModel::EVI:
    return d;  // Diagonal elements (varying across clusters)
  case CovarianceModel::VVI:
    return d;  // Full diagonal per cluster
  case CovarianceModel::FULL:
  default:
    return d * (d + 1) / 2;  // Full covariance matrix
  }
}

// Construct covariance matrix from parameters based on model
arma::mat MVNormalMixingDistribution::constructCovarianceMatrix(
    const arma::vec& params, int d) const {

  arma::mat sigma(d, d, arma::fill::zeros);

  switch (covModel) {
  case CovarianceModel::E:
  case CovarianceModel::V:
    // For univariate case, return scalar as 1x1 matrix
    sigma(0, 0) = params(0);
    break;

  case CovarianceModel::EII:
  case CovarianceModel::VII:
    // Spherical: sigma = lambda * I
    sigma = params(0) * arma::eye(d, d);
    break;

  case CovarianceModel::EEI:
    // Diagonal, equal volume and shape
    for (int i = 0; i < d; i++) {
      sigma(i, i) = params(i);
    }
    break;

  case CovarianceModel::VEI:
    // Diagonal, varying volume, equal shape
    // params[0] = volume, params[1:d] = shape
  {
    double volume = params(0);
    arma::vec shape = params.subvec(1, d);
    shape = shape / arma::prod(shape);  // Normalize shape
    for (int i = 0; i < d; i++) {
      sigma(i, i) = volume * shape(i);
    }
  }
    break;

  case CovarianceModel::EVI:
  case CovarianceModel::VVI:
    // Diagonal matrices
    for (int i = 0; i < d; i++) {
      sigma(i, i) = params(i);
    }
    break;

  case CovarianceModel::FULL:
  default:
    // Full covariance matrix (lower triangular parameterization)
  {
    int idx = 0;
    for (int i = 0; i < d; i++) {
      for (int j = 0; j <= i; j++) {
        sigma(i, j) = params(idx);
        if (i != j) sigma(j, i) = params(idx);
        idx++;
      }
    }
  }
    break;
  }

  return sigma;
}

// Extract covariance parameters from matrix based on model
arma::vec MVNormalMixingDistribution::extractCovarianceParams(
    const arma::mat& sigma) const {

  int d = sigma.n_rows;
  int nParams = getNumCovParams(d);
  arma::vec params(nParams);

  switch (covModel) {
  case CovarianceModel::E:
  case CovarianceModel::V:
    params(0) = sigma(0, 0);
    break;

  case CovarianceModel::EII:
  case CovarianceModel::VII:
    // Extract volume (average of diagonal elements)
    params(0) = arma::trace(sigma) / d;
    break;

  case CovarianceModel::EEI:
  case CovarianceModel::EVI:
  case CovarianceModel::VVI:
    // Extract diagonal elements
    for (int i = 0; i < d; i++) {
      params(i) = sigma(i, i);
    }
    break;

  case CovarianceModel::VEI:
    // Extract volume and shape
  {
    arma::vec diag = sigma.diag();
    params(0) = arma::prod(diag);  // Volume
    arma::vec shape = diag / std::pow(params(0), 1.0/d);
    params.subvec(1, d) = shape;
  }
    break;

  case CovarianceModel::FULL:
  default:
    // Extract lower triangular elements
  {
    int idx = 0;
    for (int i = 0; i < d; i++) {
      for (int j = 0; j <= i; j++) {
        params(idx) = sigma(i, j);
        idx++;
      }
    }
  }
    break;
  }

  return params;
}

arma::vec MVNormalMixingDistribution::mvnLikelihood(const arma::mat& x,
                                                    const arma::vec& mu,
                                                    const arma::mat& sigma) const {
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec result(n);

  // Ensure sigma is symmetric before inversion
  arma::mat sigma_sym = ensureSymmetric(sigma);

  // sigma here is actually a precision matrix (inverse covariance)
  // We need to convert it to covariance for likelihood calculation
  arma::mat covariance;
  try {
    covariance = arma::inv_sympd(sigma_sym);
  } catch(...) {
    result.fill(1e-300);
    return result;
  }

  // Calculate log-determinant and inverse of covariance
  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, covariance);

  if (sign <= 0) {
    // Covariance is not positive definite
    result.fill(1e-300);
    return result;
  }

  // Use the precision matrix (sigma_sym) directly for the quadratic form
  double log_const = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * log_det_val;

  for (int i = 0; i < n; i++) {
    arma::vec x_centered = x.row(i).t() - mu;
    double quad_form = arma::as_scalar(x_centered.t() * sigma_sym * x_centered);
    result(i) = std::exp(log_const - 0.5 * quad_form);
  }

  return result;
}

Rcpp::NumericVector MVNormalMixingDistribution::likelihood(const arma::vec& x,
                                                           const Rcpp::List& theta) const {
  // Extract parameters - handle the array structure
  Rcpp::NumericVector mu_array = theta["mu"];
  Rcpp::NumericVector sig_array = theta["sig"];

  // Get dimensions
  Rcpp::IntegerVector mu_dim = mu_array.attr("dim");
  int d = mu_dim[1];

  // Extract mu
  arma::vec mu(mu_array.begin(), d);

  // Handle sigma based on covariance model
  arma::mat sig;
  if (covModel == CovarianceModel::FULL) {
    // Full covariance matrix
    sig = arma::mat(sig_array.begin(), d, d);
  } else {
    // Reconstruct covariance from parameters
    int nParams = getNumCovParams(d);
    arma::vec params(sig_array.begin(), nParams);
    arma::mat cov = constructCovarianceMatrix(params, d);
    // Convert to precision
    sig = arma::inv_sympd(cov);
  }

  // Compute likelihood
  arma::mat x_mat = arma::mat(x.memptr(), 1, x.n_elem);
  arma::vec lik_vec = mvnLikelihood(x_mat, mu, sig);

  Rcpp::NumericVector result(1);
  result[0] = lik_vec(0);
  return result;
}

Rcpp::List MVNormalMixingDistribution::posteriorParameters(const arma::mat& x) const {
  int n = x.n_rows;
  int d = x.n_cols;

  // Special case: no data
  if (n == 0) {
    Rcpp::NumericVector mu0_vec = Rcpp::wrap(mu0);
    mu0_vec.attr("dim") = R_NilValue;

    return Rcpp::List::create(
      Rcpp::Named("mu_n") = mu0_vec,
      Rcpp::Named("t_n") = ensureSymmetric(Lambda),
      Rcpp::Named("Lambda_n") = ensureSymmetric(Lambda),
      Rcpp::Named("kappa_n") = kappa0,
      Rcpp::Named("nu_n") = nu
    );
  }

  // Compute sample statistics
  arma::vec x_bar = arma::mean(x, 0).t();

  // Posterior parameters for mean (same for all models)
  double kappa_n = kappa0 + n;
  arma::vec mu_n_arma = (kappa0 * mu0 + n * x_bar) / kappa_n;
  double nu_n = nu + n;

  // Compute scatter matrix based on covariance model
  arma::mat S(d, d, arma::fill::zeros);

  switch (covModel) {
  case CovarianceModel::E:
  case CovarianceModel::V:
    // Univariate case
    if (n > 1) {
      double var = arma::as_scalar(arma::var(x));
      S(0, 0) = (n - 1) * var;
    }
    break;

  case CovarianceModel::EII:
  case CovarianceModel::VII:
    // Spherical covariance
    if (n > 1) {
      arma::mat centered = x.each_row() - x_bar.t();
      double trace_S = arma::accu(centered % centered) / (n - 1);
      S = (trace_S / d) * arma::eye(d, d);
    }
    break;

  case CovarianceModel::EEI:
  case CovarianceModel::VEI:
  case CovarianceModel::EVI:
  case CovarianceModel::VVI:
    // Diagonal covariance
    if (n > 1) {
      arma::vec diag_var = arma::var(x, 0, 0).t();
      S = arma::diagmat(diag_var) * (n - 1);
    }
    break;

  case CovarianceModel::FULL:
  default:
    // Full covariance
    if (n > 1) {
      S = (n - 1) * arma::cov(x);
      S = ensureSymmetric(S);
    }
    break;
  }

  // Update Lambda (called t_n in R code)
  arma::vec diff = x_bar - mu0;
  arma::mat t_n = Lambda + S + (kappa0 * n / kappa_n) * (diff * diff.t());
  t_n = ensureSymmetric(t_n);

  // Convert arma::vec to plain Rcpp::NumericVector
  Rcpp::NumericVector mu_n_vec = Rcpp::wrap(mu_n_arma);
  mu_n_vec.attr("dim") = R_NilValue;

  return Rcpp::List::create(
    Rcpp::Named("mu_n") = mu_n_vec,
    Rcpp::Named("t_n") = t_n,
    Rcpp::Named("Lambda_n") = t_n,
    Rcpp::Named("kappa_n") = kappa_n,
    Rcpp::Named("nu_n") = nu_n
  );
}

Rcpp::List MVNormalMixingDistribution::priorDraw(int n) const {
  int d = mu0.n_elem;

  // Arrays to store results
  Rcpp::NumericVector mu_arr = Rcpp::NumericVector(Rcpp::Dimension(1, d, n));
  Rcpp::NumericVector sig_arr;

  // Determine storage size for covariance parameters
  if (covModel == CovarianceModel::FULL) {
    sig_arr = Rcpp::NumericVector(Rcpp::Dimension(d, d, n));
  } else {
    int nCovParams = getNumCovParams(d);
    sig_arr = Rcpp::NumericVector(Rcpp::Dimension(nCovParams, n));
  }

  // Ensure Lambda is symmetric
  arma::mat Lambda_sym = ensureSymmetric(Lambda);

  for (int i = 0; i < n; i++) {
    // Draw precision from Wishart
    arma::mat prec_draw = arma::wishrnd(Lambda_sym, nu);

    // Ensure the drawn precision matrix is symmetric
    prec_draw = ensureSymmetric(prec_draw);

    // Draw mu from Multivariate Normal given precision
    arma::mat cov_mu = arma::inv_sympd(prec_draw / kappa0);
    arma::vec mu_draw = arma::mvnrnd(mu0, cov_mu);

    // Store mu
    for (int j = 0; j < d; j++) {
      mu_arr[j + i * d] = mu_draw(j);
    }

    // Store covariance parameters based on model
    if (covModel == CovarianceModel::FULL) {
      // Store full precision matrix
      for (int j = 0; j < d; j++) {
        for (int k = 0; k < d; k++) {
          sig_arr[j + k * d + i * d * d] = prec_draw(j, k);
        }
      }
    } else {
      // Convert to covariance and extract model-specific parameters
      arma::mat cov_draw = arma::inv_sympd(prec_draw);
      arma::vec params = extractCovarianceParams(cov_draw);
      for (int j = 0; j < params.n_elem; j++) {
        sig_arr[j + i * params.n_elem] = params(j);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sig") = sig_arr
  );
}

Rcpp::List MVNormalMixingDistribution::posteriorDraw(const arma::mat& x, int n) const {
  // Get posterior parameters
  Rcpp::List post_params = posteriorParameters(x);

  arma::vec mu_n = Rcpp::as<arma::vec>(post_params["mu_n"]);
  arma::mat t_n = Rcpp::as<arma::mat>(post_params["t_n"]);
  double kappa_n = Rcpp::as<double>(post_params["kappa_n"]);
  double nu_n = Rcpp::as<double>(post_params["nu_n"]);

  int d = mu_n.n_elem;

  // Arrays to store results
  Rcpp::NumericVector mu_arr = Rcpp::NumericVector(Rcpp::Dimension(1, d, n));
  Rcpp::NumericVector sig_arr;

  // Determine storage size for covariance parameters
  if (covModel == CovarianceModel::FULL) {
    sig_arr = Rcpp::NumericVector(Rcpp::Dimension(d, d, n));
  } else {
    int nCovParams = getNumCovParams(d);
    sig_arr = Rcpp::NumericVector(Rcpp::Dimension(nCovParams, n));
  }

  // Ensure t_n is symmetric
  arma::mat t_n_sym = ensureSymmetric(t_n);

  for (int i = 0; i < n; i++) {
    // Draw precision from Wishart
    arma::mat prec_draw = arma::wishrnd(t_n_sym, nu_n);

    // Ensure the drawn precision matrix is symmetric
    prec_draw = ensureSymmetric(prec_draw);

    // Draw mu from Multivariate Normal given precision
    arma::mat cov_mu = arma::inv_sympd(ensureSymmetric(prec_draw / kappa_n));
    arma::vec mu_draw = arma::mvnrnd(mu_n, cov_mu);

    // Store mu
    for (int j = 0; j < d; j++) {
      mu_arr[j + i * d] = mu_draw(j);
    }

    // Store covariance parameters based on model
    if (covModel == CovarianceModel::FULL) {
      // Store full precision matrix
      for (int j = 0; j < d; j++) {
        for (int k = 0; k < d; k++) {
          sig_arr[j + k * d + i * d * d] = prec_draw(j, k);
        }
      }
    } else {
      // Convert to covariance and extract model-specific parameters
      arma::mat cov_draw = arma::inv_sympd(prec_draw);
      arma::vec params = extractCovarianceParams(cov_draw);
      for (int j = 0; j < params.n_elem; j++) {
        sig_arr[j + i * params.n_elem] = params(j);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_arr,
    Rcpp::Named("sig") = sig_arr
  );
}

Rcpp::NumericVector MVNormalMixingDistribution::predictive(const arma::mat& x) const {
  int n = x.n_rows;
  int d = x.n_cols;
  Rcpp::NumericVector result(n);

  double pi_const = std::pow(M_PI, -0.5 * d);

  for (int i = 0; i < n; i++) {
    arma::mat x_i = x.row(i);
    Rcpp::List post_params = posteriorParameters(x_i);

    arma::vec mu_n = Rcpp::as<arma::vec>(post_params["mu_n"]);
    arma::mat t_n = Rcpp::as<arma::mat>(post_params["t_n"]);
    double kappa_n = Rcpp::as<double>(post_params["kappa_n"]);
    double nu_n = Rcpp::as<double>(post_params["nu_n"]);

    // Calculate determinants
    double log_det_Lambda, log_det_t_n;
    double sign_Lambda, sign_t_n;
    arma::log_det(log_det_Lambda, sign_Lambda, Lambda);
    arma::log_det(log_det_t_n, sign_t_n, t_n);

    // Handle potential numerical issues
    if (sign_Lambda <= 0 || sign_t_n <= 0) {
      result[i] = 1e-300;
      continue;
    }

    double ratio_det = std::exp((nu / 2.0) * (log_det_Lambda - log_det_t_n));
    double ratio_kappa = std::pow(kappa0 / kappa_n, d / 2.0);

    // Compute multivariate gamma ratio
    double log_gamma_ratio = 0.0;
    for (int j = 1; j <= d; j++) {
      log_gamma_ratio += lgamma((nu_n + 1.0 - j) / 2.0) - lgamma((nu + 1.0 - j) / 2.0);
    }
    double gamma_ratio = std::exp(log_gamma_ratio);

    result[i] = pi_const * ratio_kappa * ratio_det * gamma_ratio;
  }

  return result;
}

// Static methods
Rcpp::List MVNormalMixingDistribution::priorDrawStatic(const Rcpp::List& priorParams, int n) {
  MVNormalMixingDistribution md(priorParams);
  return md.priorDraw(n);
}

Rcpp::List MVNormalMixingDistribution::posteriorDrawStatic(const Rcpp::List& priorParams,
                                                           const arma::mat& x, int n) {
  MVNormalMixingDistribution md(priorParams);
  return md.posteriorDraw(x, n);
}

// ConjugateMVNormalDP implementation
ConjugateMVNormalDP::ConjugateMVNormalDP() : mixingDistribution(nullptr), numberClusters(0) {
  // Constructor
}

ConjugateMVNormalDP::~ConjugateMVNormalDP() {
  if (mixingDistribution) {
    delete mixingDistribution;
  }
}

void ConjugateMVNormalDP::initialize(const Rcpp::List& dpObj) {
  // Extract data
  data = Rcpp::as<arma::mat>(dpObj["data"]);

  // Extract cluster labels (already 0-indexed from R wrapper)
  Rcpp::IntegerVector labels = dpObj["clusterLabels"];
  clusterLabels = arma::vec(labels.begin(), labels.size());

  // Initialize mixing distribution
  Rcpp::List mdObj = dpObj["mixingDistribution"];
  Rcpp::List priorParams = mdObj["priorParameters"];
  mixingDistribution = new MVNormalMixingDistribution(priorParams);

  // Extract cluster parameters if they exist
  if (dpObj.containsElementNamed("clusterParameters")) {
    clusterParameters = dpObj["clusterParameters"];
  }

  // Count clusters
  numberClusters = arma::max(clusterLabels) + 1;
}

Rcpp::List ConjugateMVNormalDP::updateClusterComponents() {
  int n = data.n_rows;
  double alpha = 1.0; // Default alpha, should be extracted from dpObj

  // Extract predictive array if it exists
  Rcpp::NumericVector predictiveArray(n);

  for (int i = 0; i < n; i++) {
    // Current data point
    arma::mat x_i = data.row(i);

    // Calculate predictive probabilities for existing clusters
    arma::vec probs(numberClusters + 1);

    for (int k = 0; k < numberClusters; k++) {
      // Get data points in cluster k (excluding current point)
      arma::uvec cluster_idx = arma::find((clusterLabels == k) &&
        (arma::linspace<arma::vec>(0, n-1, n) != i));

      if (cluster_idx.n_elem > 0) {
        arma::mat cluster_data = data.rows(cluster_idx);
        Rcpp::NumericVector pred = mixingDistribution->predictive(x_i);
        probs(k) = cluster_idx.n_elem * pred(0);
      } else {
        probs(k) = 0;
      }
    }

    // Probability of new cluster
    Rcpp::NumericVector pred_new = mixingDistribution->predictive(x_i);
    probs(numberClusters) = alpha * pred_new(0);

    // Normalize
    probs = probs / arma::sum(probs);

    // Sample new cluster assignment
    double u = R::runif(0, 1);
    double cumsum = 0;
    int new_label = 0;

    for (int k = 0; k <= numberClusters; k++) {
      cumsum += probs(k);
      if (u <= cumsum) {
        new_label = k;
        break;
      }
    }

    // Update cluster assignment
    clusterLabels(i) = new_label;

    // Update number of clusters if needed
    if (new_label == numberClusters) {
      numberClusters++;
    }

    // Store predictive probability
    predictiveArray[i] = pred_new(0);
  }

  // Clean up empty clusters
  arma::vec unique_labels = arma::unique(clusterLabels);
  numberClusters = unique_labels.n_elem;

  // Relabel clusters to be consecutive
  for (int k = 0; k < numberClusters; k++) {
    arma::uvec idx = arma::find(clusterLabels == unique_labels(k));
    clusterLabels.elem(idx).fill(k);
  }

  // Calculate points per cluster
  Rcpp::IntegerVector pointsPerCluster(numberClusters);
  for (int k = 0; k < numberClusters; k++) {
    pointsPerCluster[k] = arma::sum(clusterLabels == k);
  }

  return Rcpp::List::create(
    Rcpp::Named("clusterLabels") = Rcpp::IntegerVector(clusterLabels.begin(), clusterLabels.end()),
    Rcpp::Named("pointsPerCluster") = pointsPerCluster,
    Rcpp::Named("numberClusters") = numberClusters,
    Rcpp::Named("clusterParameters") = clusterParameters
  );
}

Rcpp::List ConjugateMVNormalDP::updateClusterParameters() {
  // Initialize new cluster parameters
  Rcpp::List newParams;

  int d = data.n_cols;
  int totalParams = numberClusters;

  // Determine parameter dimensions based on covariance model
  int nCovParams = mixingDistribution->getNumCovParams(d);

  Rcpp::NumericVector mu_all, sig_all;
  if (mixingDistribution->getCovarianceModel() == CovarianceModel::FULL) {
    mu_all = Rcpp::NumericVector(Rcpp::Dimension(1, d, totalParams));
    sig_all = Rcpp::NumericVector(Rcpp::Dimension(d, d, totalParams));
  } else {
    mu_all = Rcpp::NumericVector(Rcpp::Dimension(1, d, totalParams));
    sig_all = Rcpp::NumericVector(Rcpp::Dimension(nCovParams, totalParams));
  }

  for (int k = 0; k < numberClusters; k++) {
    // Get data points in cluster k
    arma::uvec cluster_idx = arma::find(clusterLabels == k);

    if (cluster_idx.n_elem > 0) {
      arma::mat cluster_data = data.rows(cluster_idx);

      // Draw parameters from posterior
      Rcpp::List params = mixingDistribution->posteriorDraw(cluster_data, 1);

      // Extract mu and sig
      Rcpp::NumericVector mu_k = params["mu"];
      Rcpp::NumericVector sig_k = params["sig"];

      // Store in arrays
      for (int j = 0; j < d; j++) {
        mu_all[j + k * d] = mu_k[j];
      }

      if (mixingDistribution->getCovarianceModel() == CovarianceModel::FULL) {
        for (int i = 0; i < d * d; i++) {
          sig_all[i + k * d * d] = sig_k[i];
        }
      } else {
        for (int i = 0; i < nCovParams; i++) {
          sig_all[i + k * nCovParams] = sig_k[i];
        }
      }
    } else {
      // Draw from prior for empty clusters
      Rcpp::List params = mixingDistribution->priorDraw(1);

      Rcpp::NumericVector mu_k = params["mu"];
      Rcpp::NumericVector sig_k = params["sig"];

      for (int j = 0; j < d; j++) {
        mu_all[j + k * d] = mu_k[j];
      }

      if (mixingDistribution->getCovarianceModel() == CovarianceModel::FULL) {
        for (int i = 0; i < d * d; i++) {
          sig_all[i + k * d * d] = sig_k[i];
        }
      } else {
        for (int i = 0; i < nCovParams; i++) {
          sig_all[i + k * nCovParams] = sig_k[i];
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_all,
    Rcpp::Named("sig") = sig_all
  );
}

// Export functions
Rcpp::List conjugate_mvnormal_cluster_component_update_cpp(const Rcpp::List& dpObj) {
  ConjugateMVNormalDP dp;
  dp.initialize(dpObj);
  return dp.updateClusterComponents();
}

Rcpp::List conjugate_mvnormal_cluster_parameter_update_cpp(const Rcpp::List& dpObj) {
  ConjugateMVNormalDP dp;
  dp.initialize(dpObj);
  return dp.updateClusterParameters();
}

} // namespace dp
