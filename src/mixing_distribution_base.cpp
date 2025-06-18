// src/mixing_distribution_base.cpp
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/gaussian_mixing.h"
#include "../inst/include/beta_mixing.h"  // ADD THIS LINE
#include <RcppArmadillo.h>

namespace dirichletprocess {

std::unique_ptr<MixingDistribution> MixingDistribution::create(
    const std::string& type,
    const Rcpp::List& params) {

  if (type == "gaussian") {
    double mu0 = Rcpp::as<double>(params["mu0"]);
    double kappa0 = Rcpp::as<double>(params["kappa0"]);
    double alpha0 = Rcpp::as<double>(params["alpha0"]);
    double beta0 = Rcpp::as<double>(params["beta0"]);
    return std::unique_ptr<MixingDistribution>(
      new GaussianMixing(mu0, kappa0, alpha0, beta0));
  } else if (type == "beta") {
    double alpha0 = Rcpp::as<double>(params["alpha0"]);
    double beta0 = Rcpp::as<double>(params["beta0"]);
    double maxT = params.containsElementNamed("maxT") ?
                  Rcpp::as<double>(params["maxT"]) : 1.0;

    return std::unique_ptr<MixingDistribution>(
      new BetaMixing(alpha0, beta0, maxT));
    }

  Rcpp::stop("Unknown mixing distribution type: " + type);
}

} // namespace dirichletprocess
