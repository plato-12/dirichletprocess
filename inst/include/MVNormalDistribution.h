// inst/include/MVNormalDistribution.h
#ifndef MVNORMAL_DISTRIBUTION_H
#define MVNORMAL_DISTRIBUTION_H

#include "DirichletProcessBase.h"
#include "mixing_distribution_base.h"
#include <RcppArmadillo.h>
#include <string>

namespace dp {

// Enum for covariance model types
enum class CovarianceModel {
  FULL,    // Full covariance matrix (default)
  E,       // Equal variance (univariate)
  V,       // Variable variance (univariate)
  EII,     // Spherical, equal volume
  VII,     // Spherical, unequal volume
  EEI,     // Diagonal, equal volume and shape
  VEI,     // Diagonal, varying volume, equal shape
  EVI,     // Diagonal, equal volume, varying shape
  VVI      // Diagonal, varying volume and shape
};

// Helper function to ensure matrix symmetry
inline arma::mat ensureSymmetric(const arma::mat& M) {
  return 0.5 * (M + M.t());
}

class MVNormalMixingDistribution : public MixingDistribution {
private:
  // Prior parameters
  arma::vec mu0;
  double kappa0;
  arma::mat Lambda;
  double nu;

  // Covariance model
  CovarianceModel covModel;

  // Helper functions for different covariance structures
  arma::mat constructCovarianceMatrix(const arma::vec& params, int d) const;
  arma::vec extractCovarianceParams(const arma::mat& sigma) const;
  // Remove getNumCovParams from here

public:
  MVNormalMixingDistribution(const Rcpp::List& priorParams);
  ~MVNormalMixingDistribution();

  // Core functions
  Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const override;
  Rcpp::List posteriorParameters(const arma::mat& x) const;
  Rcpp::List priorDraw(int n) const override;
  Rcpp::List posteriorDraw(const arma::mat& x, int n) const override;
  Rcpp::NumericVector predictive(const arma::mat& x) const;

  // MVNormal specific likelihood
  arma::vec mvnLikelihood(const arma::mat& x, const arma::vec& mu,
                          const arma::mat& sigma) const;

  // Static methods for R interface
  static Rcpp::List priorDrawStatic(const Rcpp::List& priorParams, int n);
  static Rcpp::List posteriorDrawStatic(const Rcpp::List& priorParams,
                                        const arma::mat& x, int n);

  // Get covariance model
  CovarianceModel getCovarianceModel() const { return covModel; }

  // Add this method as public
  int getNumCovParams(int d) const;
};

// Conjugate MVNormal Dirichlet Process
class ConjugateMVNormalDP {
private:
  MVNormalMixingDistribution* mixingDistribution;
  arma::mat data;
  arma::uvec clusterLabels;
  int numberClusters;
  Rcpp::List clusterParameters;
  arma::uvec pointsPerCluster;
  double alpha;
  Rcpp::List alphaPriorParameters;
  int n;
  Rcpp::NumericVector predictiveArray;

  // Internal methods
  void initialisePredictive();
  void clusterComponentUpdate();
  void clusterParameterUpdate();
  void updateAlpha();
  Rcpp::List clusterLabelChange(int i, int newLabel, int currentLabel);

public:
  ConjugateMVNormalDP();
  ~ConjugateMVNormalDP();

  void initialize(const Rcpp::List& dpObj);
  Rcpp::List updateClusterComponents();
  Rcpp::List updateClusterParameters();
};

// Export functions
Rcpp::List conjugate_mvnormal_cluster_component_update_cpp(const Rcpp::List& dpObj);
Rcpp::List conjugate_mvnormal_cluster_parameter_update_cpp(const Rcpp::List& dpObj);

} // namespace dp

#endif
