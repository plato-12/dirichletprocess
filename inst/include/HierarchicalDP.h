// inst/include/HierarchicalDP.h
#ifndef HIERARCHICAL_DP_H
#define HIERARCHICAL_DP_H

#include "DirichletProcess.h"
#include <vector>

namespace dp {

class HierarchicalDP : public DirichletProcess {
public:
  HierarchicalDP();
  virtual ~HierarchicalDP();

  // Hierarchical DP specific attributes
  std::vector<DirichletProcess*> indDP;  // Individual DPs
  Rcpp::List globalParameters;           // Global parameters
  arma::vec globalStick;                 // Global stick-breaking weights
  double gamma;                          // Concentration parameter for the top level DP
  Rcpp::NumericVector gammaPriors;       // Prior parameters for gamma

  // Implementation of core MCMC methods for hierarchical case
  void clusterComponentUpdate() override;
  void clusterParameterUpdate() override;
  void updateAlpha() override;

  // Hierarchical specific methods
  void globalParameterUpdate();
  void updateG0();
  void updateGamma();

  // Conversion methods
  Rcpp::List toR() const override;
  static HierarchicalDP* fromR(const Rcpp::List& rObj);
};

// Specialized classes for specific distribution types
class HierarchicalBetaDP : public HierarchicalDP {
public:
  HierarchicalBetaDP();
  virtual ~HierarchicalBetaDP();

  // Additional methods specific to Beta hierarchical DP
};

class HierarchicalMVNormal2DP : public HierarchicalDP {
public:
  HierarchicalMVNormal2DP();
  virtual ~HierarchicalMVNormal2DP();

  // Additional methods specific to MVNormal2 hierarchical DP
};

} // namespace dp

#endif
