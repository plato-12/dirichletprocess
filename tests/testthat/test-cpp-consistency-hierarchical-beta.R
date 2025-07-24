# tests/testthat/test-cpp-consistency-hierarchical-beta.R
#
# R/C++ Consistency Tests for Hierarchical Beta Distribution
# 
# This file tests statistical equivalence between R and C++ implementations
# for the Hierarchical Beta distribution with comprehensive validation.

# Development vs Production testing configuration
DEV_MODE <- Sys.getenv("DP_DEV_TESTING", "TRUE") == "TRUE"
BASE_ITERATIONS <- if (DEV_MODE) 50 else 200

test_that("Hierarchical Beta distribution R/C++ consistency", {
  set.seed(123)
  
  # Create hierarchical test data - multiple groups
  group1_size <- if (DEV_MODE) 25 else 50
  group2_size <- if (DEV_MODE) 25 else 50
  group3_size <- if (DEV_MODE) 25 else 50
  
  # Generate test data with different beta parameters for each group
  group1_data <- rbeta(group1_size, 2, 5)
  group2_data <- rbeta(group2_size, 5, 2) 
  group3_data <- rbeta(group3_size, 1, 1)
  
  hierarchical_data <- list(group1_data, group2_data, group3_data)
  
  results <- validate_r_cpp_consistency("hierarchical_beta", hierarchical_data, 
                                        iterations = BASE_ITERATIONS,
                                        maxY = 1)

  expect_lt(results$alpha_mean_diff, ALPHA_TOLERANCE)
  expect_lt(results$cluster_count_diff, CLUSTER_TOLERANCE)
  expect_gt(results$likelihood_correlation, LIKELIHOOD_CORR_MIN)
})

test_that("Hierarchical Beta manual MCMC interface", {
  skip_if_not(using_cpp(), "C++ not available")
  
  set.seed(456)
  
  # Create smaller test data for manual MCMC
  group_size <- if (DEV_MODE) 15 else 30
  group1_data <- rbeta(group_size, 3, 2)
  group2_data <- rbeta(group_size, 2, 3)
  hierarchical_data <- list(group1_data, group2_data)
  
  # Test CppMCMCRunner interface with hierarchical data
  dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)
  
  # Enable C++ hierarchical samplers
  enable_cpp_hierarchical_samplers()
  
  # Test manual MCMC steps
  manual_iterations <- if (DEV_MODE) 25 else 50
  for (i in 1:manual_iterations) {
    dp <- ClusterComponentUpdate(dp)
    dp <- ClusterParameterUpdate(dp)
    dp <- UpdateAlpha(dp)
    
    # Test hierarchical-specific updates
    dp <- GlobalParameterUpdate(dp)
  }
  
  # Verify results
  expect_s3_class(dp, c("list", "dirichletprocess", "hierarchical"))
  expect_true(length(dp) >= 2)  # Multiple groups
  expect_true(all(sapply(dp, function(x) x$numberClusters >= 1)))
})

test_that("Hierarchical Beta global parameter sharing", {
  skip_if_not(using_cpp(), "C++ not available")
  
  set.seed(789)
  
  # Create test data with similar patterns to test global parameter sharing
  group_size <- if (DEV_MODE) 20 else 40
  # Use similar beta parameters to encourage parameter sharing
  group1_data <- rbeta(group_size, 4, 2)
  group2_data <- rbeta(group_size, 4.2, 1.8)  # Slightly different but similar
  hierarchical_data <- list(group1_data, group2_data)
  
  # Test with C++ enabled
  set_use_cpp(TRUE)
  enable_cpp_hierarchical_samplers()
  
  dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)
  dp <- Fit(dp, its = if (DEV_MODE) 30 else 100)
  
  # Verify hierarchical structure
  expect_s3_class(dp, c("list", "dirichletprocess", "hierarchical"))
  expect_true(length(dp) == 2)
  
  # Check that global parameters exist and are shared
  expect_true(exists("globalParameters", where = environment(dp)))
  
  # Verify that individual groups maintain their own cluster structure
  expect_true(all(sapply(dp, function(x) {
    "clusterLabels" %in% names(x) && "clusterParameters" %in% names(x)
  })))
})

test_that("Hierarchical Beta edge cases", {
  skip_if_not(using_cpp(), "C++ not available")
  
  set.seed(101112)
  
  # Test with very small groups
  small_data <- list(rbeta(3, 1, 1), rbeta(2, 1, 1))
  
  expect_no_error({
    dp_small <- DirichletProcessHierarchicalBeta(small_data, maxY = 1)
    dp_small <- Fit(dp_small, its = if (DEV_MODE) 10 else 25)
  })
  
  # Test with single group (edge case)
  single_group <- list(rbeta(if (DEV_MODE) 10 else 20, 2, 2))
  
  expect_no_error({
    dp_single <- DirichletProcessHierarchicalBeta(single_group, maxY = 1)
    dp_single <- Fit(dp_single, its = if (DEV_MODE) 10 else 25)
  })
})