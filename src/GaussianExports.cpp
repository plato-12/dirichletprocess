// src/GaussianExports.cpp
#include <RcppArmadillo.h>
#include "../inst/include/mcmc_runner.h"
#include "../inst/include/mixing_distribution_base.h"

// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(arma::mat data,
                        Rcpp::List mixing_dist_params,
                        Rcpp::List mcmc_params) {
  try {
    // Input validation
    if (data.n_rows == 0 || data.n_cols == 0) {
      Rcpp::stop("Data matrix cannot be empty");
    }

    if (data.has_nan()) {
      Rcpp::stop("Data contains NA values");
    }

    if (data.has_inf()) {
      Rcpp::stop("Data contains infinite values");
    }

    dirichletprocess::MCMCRunner runner(data, mixing_dist_params, mcmc_params);
    return runner.run();

  } catch (const std::exception& e) {
    Rcpp::stop("C++ MCMC error: " + std::string(e.what()));
  } catch (...) {
    Rcpp::stop("Unknown error in C++ MCMC");
  }
}
