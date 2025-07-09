// inst/include/mixing_distribution_base.h
#ifndef MIXING_DISTRIBUTION_BASE_H
#define MIXING_DISTRIBUTION_BASE_H

#include <RcppArmadillo.h>
#include <memory>

namespace dirichletprocess {

class MixingDistribution {
public:
  virtual ~MixingDistribution() {}

  // Log likelihood of data point given parameters
  virtual double log_likelihood(const arma::vec& data_point,
                                const arma::vec& params) const = 0;

  // Sample from posterior given cluster data
  virtual arma::vec posterior_draw(const arma::mat& cluster_data,
                                   const arma::vec& prior_params) const = 0;

  // Sample from prior
  virtual arma::vec prior_draw() const = 0;

  // Parameter dimension
  virtual int param_dim() const = 0;

  // Factory method
  static std::unique_ptr<MixingDistribution> create(const std::string& type,
                                                    const Rcpp::List& params);
};

} // namespace dirichletprocess

#endif
