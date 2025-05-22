// inst/include/DirichletProcessBase.h
#ifndef DIRICHLETPROCESS_BASE_H
#define DIRICHLETPROCESS_BASE_H

#include <RcppArmadillo.h>

namespace dp {

// Forward declarations
class MixingDistribution;
class DirichletProcess;

// Base MixingDistribution class
class MixingDistribution {
public:
  MixingDistribution();
  virtual ~MixingDistribution();

  // Core properties
  std::string distribution;
  bool conjugate;
  Rcpp::RObject priorParameters;
  Rcpp::RObject mhStepSize;
  Rcpp::RObject hyperPriorParameters;

  // Virtual methods that should be implemented by derived classes (with default implementations)
  virtual Rcpp::NumericVector likelihood(const arma::vec& x, const Rcpp::List& theta) const {
    Rcpp::stop("Base MixingDistribution::likelihood() called - must be overridden");
    return Rcpp::NumericVector();
  }
  virtual Rcpp::List priorDraw(int n) const {
    Rcpp::stop("Base MixingDistribution::priorDraw() called - must be overridden");
    return Rcpp::List();
  }
  virtual Rcpp::List posteriorDraw(const arma::mat& x, int n = 1) const {
    Rcpp::stop("Base MixingDistribution::posteriorDraw() called - must be overridden");
    return Rcpp::List();
  }

  // Conversion methods
  virtual Rcpp::List toR() const;
};

// Base DirichletProcess class
class DirichletProcess {
public:
  DirichletProcess();
  virtual ~DirichletProcess();

  // Common properties
  arma::mat data;
  int n;
  double alpha;
  Rcpp::RObject alphaPriorParameters;

  // Virtual methods for MCMC algorithms (with default implementations)
  virtual void clusterComponentUpdate() {
    Rcpp::warning("Base DirichletProcess::clusterComponentUpdate() called - should be overridden");
  }
  virtual void clusterParameterUpdate() {
    Rcpp::warning("Base DirichletProcess::clusterParameterUpdate() called - should be overridden");
  }
  virtual void updateAlpha() {
    Rcpp::warning("Base DirichletProcess::updateAlpha() called - should be overridden");
  }

  // Conversion methods
  virtual Rcpp::List toR() const;
  static DirichletProcess* fromR(const Rcpp::List& rObj);
};

} // namespace dp

#endif
