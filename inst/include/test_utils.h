// File: inst/include/test_utils.h

#ifndef DIRICHLETPROCESS_TEST_UTILS_H
#define DIRICHLETPROCESS_TEST_UTILS_H

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"       // Include Catch2 header
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
#include <random>

namespace DPTest {

// Constants for testing
const double EPSILON = 1e-10;  // Tolerance for floating point comparisons

// Utility function to generate random data for testing
inline Rcpp::NumericVector generateRandomData(int n, double mean = 0.0, double sd = 1.0) {
  Rcpp::NumericVector result(n);
  std::mt19937 gen(42);  // Fixed seed for reproducibility
  std::normal_distribution<double> dist(mean, sd);

  for (int i = 0; i < n; i++) {
    result[i] = dist(gen);
  }
  return result;
}

// Create a simplified mixing distribution object for testing
inline Rcpp::List createTestMixingDistribution(const std::string& distributionType,
                                               const Rcpp::NumericVector& priorParams) {
  Rcpp::List mdObj;
  mdObj["distribution"] = distributionType;
  mdObj["priorParameters"] = priorParams;
  mdObj["conjugate"] = true;  // Default to conjugate for testing

  return mdObj;
}

// Create parameter objects for different distributions
inline Rcpp::List createNormalParameters(double mean, double sd) {
  Rcpp::List theta = Rcpp::List::create(mean, sd);
  return theta;
}

inline Rcpp::List createBetaParameters(double alpha, double beta, double maxT = 1.0) {
  Rcpp::List theta = Rcpp::List::create(
    Rcpp::Named("mu") = Rcpp::NumericVector::create(alpha/(alpha+beta) * maxT),
    Rcpp::Named("nu") = Rcpp::NumericVector::create((alpha+beta))
  );
  return theta;
}

// Utility for comparing vectors with tolerance
inline bool vectorsEqual(const Rcpp::NumericVector& a,
                         const Rcpp::NumericVector& b,
                         double epsilon = EPSILON) {
  if (a.size() != b.size()) return false;

  for (int i = 0; i < a.size(); i++) {
    if (std::abs(a[i] - b[i]) > epsilon) return false;
  }
  return true;
}

// Calculate R-equivalent normal likelihood (for validation)
inline Rcpp::NumericVector rnorm_density(const Rcpp::NumericVector& x,
                                         double mean,
                                         double sd) {
  int n = x.size();
  Rcpp::NumericVector result(n);

  double normalizer = 1.0 / (sd * std::sqrt(2.0 * M_PI));

  for (int i = 0; i < n; i++) {
    double z = (x[i] - mean) / sd;
    result[i] = normalizer * std::exp(-0.5 * z * z);
  }

  return result;
}
}

#endif // DIRICHLETPROCESS_TEST_UTILS_H
