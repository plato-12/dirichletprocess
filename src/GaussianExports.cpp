// src/GaussianExports.cpp
#include <RcppArmadillo.h>
#include "../inst/include/mcmc_runner.h"

// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(arma::mat data,
                        Rcpp::List mixing_dist_params,
                        Rcpp::List mcmc_params) {
  dirichletprocess::MCMCRunner runner(data, mixing_dist_params, mcmc_params);
  return runner.run();
}
