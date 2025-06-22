// src/mixing_distribution_base.cpp
#include "../inst/include/mixing_distribution_base.h"
#include "../inst/include/gaussian_mixing.h"
#include "../inst/include/beta_mixing.h"
#include "../inst/include/mvnormal_mixing.h"
#include "../inst/include/weibull_mixing.h"
#include "../inst/include/exponential_mixing.h"
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
  } else if (type == "mvnormal") {
    arma::vec mu0 = Rcpp::as<arma::vec>(params["mu0"]);
    double kappa0 = Rcpp::as<double>(params["kappa0"]);
    arma::mat Lambda = Rcpp::as<arma::mat>(params["Lambda"]);
    double nu = Rcpp::as<double>(params["nu"]);

    return std::unique_ptr<MixingDistribution>(
      new MVNormalMixing(mu0, kappa0, Lambda, nu));
  } else if (type == "weibull") {
    double phi = Rcpp::as<double>(params["phi"]);
    double alpha0 = Rcpp::as<double>(params["alpha0"]);
    double beta0 = Rcpp::as<double>(params["beta0"]);

    // Optional hyperprior parameters
    double hyper_a1 = params.containsElementNamed("hyper_a1") ?
    Rcpp::as<double>(params["hyper_a1"]) : 6.0;
    double hyper_a2 = params.containsElementNamed("hyper_a2") ?
    Rcpp::as<double>(params["hyper_a2"]) : 2.0;
    double hyper_b1 = params.containsElementNamed("hyper_b1") ?
    Rcpp::as<double>(params["hyper_b1"]) : 1.0;
    double hyper_b2 = params.containsElementNamed("hyper_b2") ?
    Rcpp::as<double>(params["hyper_b2"]) : 0.5;

    // MH parameters
    double mh_step_alpha = params.containsElementNamed("mh_step_alpha") ?
    Rcpp::as<double>(params["mh_step_alpha"]) : 0.1;
    int mh_draws = params.containsElementNamed("mh_draws") ?
    Rcpp::as<int>(params["mh_draws"]) : 100;

    return std::unique_ptr<MixingDistribution>(
      new WeibullMixing(phi, alpha0, beta0, hyper_a1, hyper_a2,
                        hyper_b1, hyper_b2, mh_step_alpha, mh_draws));
  } else if (type == "exponential") {
    // Extract prior parameters
    double alpha0 = params["alpha0"];
    double beta0 = params["beta0"];

    return std::make_unique<ExponentialMixing>(alpha0, beta0);
  }


  Rcpp::stop("Unknown mixing distribution type: " + type);
}

} // namespace dirichletprocess
