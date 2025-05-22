// inst/include/MVNormalDistribution.h
#ifndef MVNORMAL_DISTRIBUTION_H
#define MVNORMAL_DISTRIBUTION_H

#include "DirichletProcess.h"

namespace dp {

class MVNormalMixingDistribution : public MixingDistribution {
public:
  MVNormalMixingDistribution(const Rcpp::List& priorParams);
  virtual ~MVNormalMixingDistribution();

  // Prior parameters specific to MVN
  arma::rowvec mu0;
  arma::mat Lambda;
  double kappa0;
  double nu;

  // Implement required virtual methods
  Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const override;
  arma::vec mvnLikelihood(const arma::mat& x, const arma::rowvec& mu, const arma::mat& sigma) const;
  Rcpp::List priorDraw(int n) const override;
  Rcpp::List posteriorDraw(const arma::mat& x, int n = 1) const override;

  // Specific methods for MVN distribution
  Rcpp::List posteriorParameters(const arma::mat& x) const;
  Rcpp::NumericVector predictive(const arma::mat& x) const;
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

  // Implementation of core MCMC methods
  void clusterComponentUpdate() override;
  void clusterParameterUpdate() override;
  void updateAlpha() override;

  // Additional methods
  Rcpp::List clusterLabelChange(int i, int newLabel, int currentLabel);
};

} // namespace dp

#endif
