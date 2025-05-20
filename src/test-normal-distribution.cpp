// File: src/test-normal-distribution.cpp

#include "../inst/include/test_utils.h"
#include "likelihood_functions.h"

// Test fixture for normal distribution tests
class NormalDistributionTest {
public:
  NormalDistributionTest() :
  mdObj(DPTest::createTestMixingDistribution("normal", Rcpp::NumericVector::create(0, 1, 1, 1))),
  smallData(DPTest::generateRandomData(10)),
  largeData(DPTest::generateRandomData(1000)) {
  }

  Rcpp::List mdObj;
  Rcpp::NumericVector smallData;
  Rcpp::NumericVector largeData;
};

TEST_CASE_METHOD(NormalDistributionTest, "Normal likelihood calculation", "[normal]") {
  // Test with standard normal parameters
  Rcpp::List theta = DPTest::createNormalParameters(0.0, 1.0);

  // Call the C++ implementation
  Rcpp::NumericVector result = likelihood_normal_cpp(mdObj, smallData, theta);

  // Calculate reference values
  Rcpp::NumericVector expected = DPTest::rnorm_density(smallData, 0.0, 1.0);

  // Check results
  REQUIRE(DPTest::vectorsEqual(result, expected));

  // Test with different parameters
  Rcpp::List theta2 = DPTest::createNormalParameters(2.5, 0.5);
  Rcpp::NumericVector result2 = likelihood_normal_cpp(mdObj, smallData, theta2);
  Rcpp::NumericVector expected2 = DPTest::rnorm_density(smallData, 2.5, 0.5);

  REQUIRE(DPTest::vectorsEqual(result2, expected2));
}

TEST_CASE_METHOD(NormalDistributionTest, "Normal likelihood with edge cases", "[normal]") {
  // Test with very small standard deviation
  Rcpp::List theta = DPTest::createNormalParameters(0.0, 1e-5);
  Rcpp::NumericVector result = likelihood_normal_cpp(mdObj, Rcpp::NumericVector::create(0.0), theta);

  // The density should be very high at the mean
  CHECK(result[0] > 1e4);

  // Test with data far from the mean
  Rcpp::List theta2 = DPTest::createNormalParameters(0.0, 1.0);
  Rcpp::NumericVector farData = Rcpp::NumericVector::create(10.0, -10.0);
  Rcpp::NumericVector result2 = likelihood_normal_cpp(mdObj, farData, theta2);

  // The density should be very small far from the mean
  CHECK(result2[0] < 1e-10);
  CHECK(result2[1] < 1e-10);
}
