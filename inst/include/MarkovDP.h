// inst/include/MarkovDP.h
#ifndef MARKOV_DP_H
#define MARKOV_DP_H

#include "DirichletProcess.h"

namespace dp {

class MarkovDP : public DirichletProcess {
public:
  MarkovDP();
  virtual ~MarkovDP();

  // Markov DP specific attributes
  arma::uvec states;            // State sequence
  Rcpp::List uniqueParams;      // Unique state parameters
  std::vector<Rcpp::List> params; // Parameters for each state
  double beta;                  // Transition concentration parameter

  // Implementation of core MCMC methods for Markov case
  void clusterComponentUpdate() override;
  void clusterParameterUpdate() override;
  void updateAlpha() override;

  // Markov specific methods
  void updateStates();
  void updateAlphaBeta();
  void paramUpdate();

  // Conversion methods
  Rcpp::List toR() const override;
  static MarkovDP* fromR(const Rcpp::List& rObj);
};

} // namespace dp

#endif
