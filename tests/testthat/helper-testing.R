# Helper functions for testing framework
# tests/testthat/helper-testing.R

# Generate appropriate test data for each distribution
generate_test_data <- function(distribution, n = 100) {
  set.seed(123)  # For reproducibility

  switch(distribution,
         "normal" = {
           # Mixture of 3 Gaussians
           c(rnorm(n/3, mean = -2, sd = 0.5),
             rnorm(n/3, mean = 0, sd = 1),
             rnorm(n/3, mean = 2, sd = 0.5))
         },
         "exponential" = {
           # Mixture of exponentials
           c(rexp(n/2, rate = 0.5),
             rexp(n/2, rate = 2))
         },
         "beta" = {
           # Mixture of beta distributions
           c(rbeta(n/2, shape1 = 2, shape2 = 5),
             rbeta(n/2, shape1 = 5, shape2 = 2))
         },
         "weibull" = {
           # Mixture of Weibull distributions
           c(rweibull(n/3, shape = 0.5, scale = 1),
             rweibull(n/3, shape = 1.5, scale = 1),
             rweibull(n/3, shape = 3, scale = 1))
         },
         "mvnormal" = {
           # Mixture of 2D Gaussians
           mu1 <- c(0, 0)
           mu2 <- c(3, 3)
           sigma <- diag(2)
           rbind(mvtnorm::rmvnorm(n/2, mu1, sigma),
                 mvtnorm::rmvnorm(n/2, mu2, sigma))
         },
         "mvnormal2" = {
           # Mixture with correlated covariance
           mu1 <- c(-2, -2)
           mu2 <- c(2, 2)
           sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
           rbind(mvtnorm::rmvnorm(n/2, mu1, sigma),
                 mvtnorm::rmvnorm(n/2, mu2, sigma))
         },
         stop("Unknown distribution: ", distribution)
  )
}

# Create appropriate DP object for each distribution
create_dp_object <- function(distribution, data, ...) {
  switch(distribution,
         "normal" = DirichletProcessGaussian(data, ...),
         "exponential" = DirichletProcessExponential(data, ...),
         "beta" = DirichletProcessBeta(data, ...),
         "weibull" = DirichletProcessWeibull(data, ...),
         "mvnormal" = DirichletProcessMvnormal(data, ...),
         "mvnormal2" = DirichletProcessMvnormal2(data, ...),
         stop("Unknown distribution: ", distribution)
  )
}

# Run consistency tests for all distributions
run_consistency_tests <- function() {
  distributions <- c("normal", "exponential", "beta", "weibull", "mvnormal", "mvnormal2")
  results <- list()

  for (dist in distributions) {
    cat("Testing", dist, "consistency...\n")
    test_data <- generate_test_data(dist, n = 100)
    results[[dist]] <- validate_r_cpp_consistency(dist, test_data, iterations = 100)
  }

  return(results)
}

# Run edge case tests
run_edge_case_tests <- function() {
  test_files <- list.files("tests/testthat", pattern = "test-cpp-edge-cases", full.names = TRUE)
  testthat::test_file(test_files)
}

# Run convergence tests
run_convergence_tests <- function() {
  test_files <- list.files("tests/testthat", pattern = "test-cpp-convergence", full.names = TRUE)
  testthat::test_file(test_files)
}
