// inst/include/BetaDistribution.h
#ifndef BETA_DISTRIBUTION_H
#define BETA_DISTRIBUTION_H

#include "DirichletProcess.h"

namespace dp {

class BetaMixingDistribution : public MixingDistribution {
public:
  BetaMixingDistribution(const Rcpp::NumericVector& priorParams);
  virtual ~BetaMixingDistribution();

  double maxT; // Upper bound for the Beta distribution

  // Implement required virtual methods
  Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const override;
  Rcpp::List priorDraw(int n) const override;
  Rcpp::List posteriorDraw(const arma::mat& x, int n = 1) const override;

  // Specific methods for Beta distribution
  Rcpp::NumericVector priorDensity(const Rcpp::List& theta) const;
  Rcpp::List mhParameterProposal(const Rcpp::List& oldParams) const;
  Rcpp::List penalisedLikelihood(const arma::mat& x) const;
  void updatePriorParameters(const Rcpp::List& clusterParameters, int n = 1);
};

class ConjugateBetaDP : public DirichletProcess {
public:
  ConjugateBetaDP();
  virtual ~ConjugateBetaDP();

  BetaMixingDistribution* mixingDistribution;

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

class NonConjugateBetaDP : public DirichletProcess {
public:
  NonConjugateBetaDP();
  virtual ~NonConjugateBetaDP();

  BetaMixingDistribution* mixingDistribution;

  // Cluster information
  arma::uvec clusterLabels;
  arma::uvec pointsPerCluster;
  int numberClusters;
  Rcpp::List clusterParameters;
  int m; // Number of auxiliary parameters

  // Implementation of core MCMC methods
  void clusterComponentUpdate() override;
  void clusterParameterUpdate() override;
  void updateAlpha() override;

  // Additional methods specific to non-conjugate Beta
  Rcpp::List clusterLabelChange(int i, int newLabel, int currentLabel, const Rcpp::List& aux);
  Rcpp::List metropolisHastings(const arma::mat& x, const Rcpp::List& startPos, int noDraws);
};

} // namespace dp

#endif
