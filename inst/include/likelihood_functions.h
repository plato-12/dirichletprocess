// likelihood_functions.h
#ifndef LIKELIHOOD_FUNCTIONS_H
#define LIKELIHOOD_FUNCTIONS_H

#include <RcppArmadillo.h>

// Function declarations
Rcpp::NumericVector likelihood_normal_cpp(
    const Rcpp::List& mdObj,
    const Rcpp::NumericVector& x,
    const Rcpp::List& theta);

#endif // LIKELIHOOD_FUNCTIONS_H
