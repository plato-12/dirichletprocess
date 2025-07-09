// inst/include/MVNormalDistribution.h
#ifndef MVNORMAL_DISTRIBUTION_H
#define MVNORMAL_DISTRIBUTION_H

#include "DirichletProcessBase.h"
#include <RcppArmadillo.h>

namespace dp {

// Helper function to ensure matrix symmetry
inline arma::mat ensureSymmetric(const arma::mat& A) {
  return 0.5 * (A + A.t());
}

class MVNormalMixingDistribution : public MixingDistribution {
public:
  MVNormalMixingDistribution(const Rcpp::List& priorParams);
  virtual ~MVNormalMixingDistribution();

  // Prior parameters specific to MVN-Wishart
  arma::vec mu0;      // Prior mean vector
  double kappa0;      // Prior precision parameter for mean
  arma::mat Lambda;   // Prior scale matrix (inverse of Psi in some texts)
  double nu;          // Prior degrees of freedom

  // Implement required virtual methods
  Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const override;
  Rcpp::List priorDraw(int n) const override;
  Rcpp::List posteriorDraw(const arma::mat& x, int n = 1) const override;

  // Helper method for multivariate normal likelihood calculation
  arma::vec mvnLikelihood(const arma::mat& x, const arma::vec& mu, const arma::mat& sigma) const;

  // Specific methods for MVN distribution
  Rcpp::List posteriorParameters(const arma::mat& x) const;
  Rcpp::NumericVector predictive(const arma::mat& x) const;

  // Static methods for direct testing
  static Rcpp::List priorDrawStatic(const Rcpp::List& priorParams, int n);
  static Rcpp::List posteriorDrawStatic(const Rcpp::List& priorParams, const arma::mat& x, int n);

private:

};

class ConjugateMVNormalDP : public DirichletProcess {
public:
  ConjugateMVNormalDP();
  virtual ~ConjugateMVNormalDP();

  MVNormalMixingDistribution* mixingDistribution;

  // Cluster information
  arma::uvec clusterLabels;
  arma::uvec pointsPerCluster;
  int numberClusters;
  Rcpp::List clusterParameters;
  arma::vec predictiveArray;

  // Implementation of core MCMC methods
  void clusterComponentUpdate() override;
  void clusterParameterUpdate() override;
  void updateAlpha() override;

  // Additional methods
  Rcpp::List clusterLabelChange(int i, int newLabel, int currentLabel);
  void initialisePredictive();
};

} // namespace dp

#endif
