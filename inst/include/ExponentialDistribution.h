// inst/include/ExponentialDistribution.h
#ifndef EXPONENTIAL_DISTRIBUTION_H
#define EXPONENTIAL_DISTRIBUTION_H

#include "DirichletProcess.h"

namespace dp {

class ExponentialMixingDistribution : public MixingDistribution {
public:
  ExponentialMixingDistribution(const Rcpp::NumericVector& priorParams);
  virtual ~ExponentialMixingDistribution();

  // Implement required virtual methods
  Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const override;
  Rcpp::List priorDraw(int n) const override;
  Rcpp::List posteriorDraw(const arma::mat& x, int n = 1) const override;

  // Exponential specific methods
  Rcpp::NumericVector predictive(const arma::vec& x) const;
  Rcpp::NumericMatrix posteriorParameters(const arma::mat& x) const;
};

class ConjugateExponentialDP : public DirichletProcess {
public:
  ConjugateExponentialDP(Rcpp::List dpObj); // Updated constructor
  virtual ~ConjugateExponentialDP();

  ExponentialMixingDistribution* mixingDistribution;

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
