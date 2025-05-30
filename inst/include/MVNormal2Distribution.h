// inst/include/MVNormal2Distribution.h
#ifndef MVNORMAL2_DISTRIBUTION_H
#define MVNORMAL2_DISTRIBUTION_H

#include <RcppArmadillo.h> // For Rcpp::List, Rcpp::NumericVector, arma::rowvec, arma::mat
#include "DirichletProcess.h" // Contains definition for MixingDistribution and DirichletProcess

namespace dp {

class MVNormal2MixingDistribution : public MixingDistribution {
public:
  MVNormal2MixingDistribution(const Rcpp::List& priorParams);
  virtual ~MVNormal2MixingDistribution();

  // ---- ADDED/MODIFIED MEMBERS ----
  std::string distribution;       // Was used in .cpp, now declared
  bool conjugate;                // Was used in .cpp, now declared
  Rcpp::List priorParameters;    // CRITICAL: Declare as Rcpp::List. This was the source of the error.
  Rcpp::NumericVector mhStepSize;// Was used in .cpp, now declared
  // ------------------------------

  // Prior parameters specific to MVNormal2 (were already correctly declared)
  arma::rowvec mu0;
  arma::mat sigma0;
  arma::mat phi0;
  double nu0;

  // Implement required virtual methods (signatures were okay but depend on Rcpp types)
  Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const override;
  Rcpp::List priorDraw(int n) const override;
  Rcpp::List posteriorDraw(const arma::mat& x, int n = 1) const override;
};

class NonConjugateMVNormal2DP : public DirichletProcess {
public:
  NonConjugateMVNormal2DP();
  virtual ~NonConjugateMVNormal2DP();

  // This member type (pointer to MVNormal2MixingDistribution) is fine.
  MVNormal2MixingDistribution* mixingDistribution;

  // Cluster information
  arma::uvec clusterLabels;
  arma::uvec pointsPerCluster;
  int numberClusters;
  Rcpp::List clusterParameters; // Uses Rcpp::List, needs <RcppArmadillo.h>
  int m; // Number of auxiliary parameters

  // Implementation of core MCMC methods
  void clusterComponentUpdate() override;
  void clusterParameterUpdate() override;
  void updateAlpha() override;

  // Additional methods
  Rcpp::List clusterLabelChange(int i, int newLabel, int currentLabel, const Rcpp::List& aux); // Uses Rcpp::List
};

} // namespace dp

#endif // MVNORMAL2_DISTRIBUTION_H
