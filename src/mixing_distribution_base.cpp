// src/mixing_distribution_base.cpp
#include <RcppArmadillo.h>
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/gaussian_mixing.h"
#include <memory>

namespace dirichletprocess {

std::unique_ptr<MixingDistribution> MixingDistribution::create(
    const std::string& type,
    const Rcpp::List& params) {

  if (type == "gaussian") {
    double mu0 = Rcpp::as<double>(params["mu0"]);
    double kappa0 = Rcpp::as<double>(params["kappa0"]);
    double alpha0 = Rcpp::as<double>(params["alpha0"]);
    double beta0 = Rcpp::as<double>(params["beta0"]);

    // C++11 compatible way to create unique_ptr
    return std::unique_ptr<MixingDistribution>(
      new GaussianMixing(mu0, kappa0, alpha0, beta0));
  }

  Rcpp::stop("Unknown mixing distribution type: " + type);
  return nullptr;
}

} // namespace dirichletprocess
