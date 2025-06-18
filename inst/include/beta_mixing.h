#ifndef BETA_MIXING_H
#define BETA_MIXING_H

#include "mixing_distribution_base.h"
#include <RcppArmadillo.h>
#include <limits>
#include <algorithm>

namespace dirichletprocess {

class BetaMixing : public MixingDistribution {
private:
  double alpha0;  // Prior shape parameter for alpha
  double beta0;   // Prior shape parameter for beta
  double maxT;    // Upper bound of Beta distribution (default 1)

  // Method of moments helper function
  arma::vec method_of_moments_estimate(const arma::vec& x) const;

public:
  BetaMixing(double alpha0, double beta0, double maxT = 1.0)
    : alpha0(alpha0), beta0(beta0), maxT(maxT) {}

  double log_likelihood(const arma::vec& data_point,
                        const arma::vec& params) const override;

  arma::vec posterior_draw(const arma::mat& cluster_data,
                           const arma::vec& prior_params) const override;

  arma::vec prior_draw() const override;

  int param_dim() const override { return 2; }  // [mu, tau] parameterization

  // Additional methods for non-conjugate MCMC
  arma::vec metropolis_hastings_step(const arma::mat& cluster_data,
                                     const arma::vec& current_params,
                                     const arma::vec& step_sizes,
                                     int n_draws = 250) const;
};

} // namespace dirichletprocess

#endif // BETA_MIXING_H
