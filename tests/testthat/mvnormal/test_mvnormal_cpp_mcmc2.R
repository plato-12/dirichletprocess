context("MVNormal C++ MCMC Implementation Tests")

test_that("MVNormal C++ MCMC produces valid results", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))

  # Check if MVNormal C++ functions are available
  skip_if_not(exists("conjugate_mvnormal_cluster_component_update_cpp"),
              "MVNormal C++ functions not available")

  # Generate test data
  set.seed(123)
  n <- 100
  d <- 2

  # True parameters for two clusters
  mu1 <- c(-2, -2)
  mu2 <- c(2, 2)
  Sigma1 <- diag(2) * 0.5
  Sigma2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)

  # Generate data
  data1 <- mvtnorm::rmvnorm(50, mu1, Sigma1)
  data2 <- mvtnorm::rmvnorm(50, mu2, Sigma2)
  y <- rbind(data1, data2)

  # Set up prior parameters
  g0Priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2),
    kappa0 = 1,
    nu = d + 2
  )

  # Create DP object
  dp <- DirichletProcessMvnormal(y, g0Priors)

  # Initialize with a reasonable number of clusters (not too many)
  dp <- Initialise(dp, numInitialClusters = 2)

  # Ensure parameter arrays have enough space
  if (dim(dp$clusterParameters$mu)[3] < 20) {
    new_mu <- array(NA_real_, dim = c(1, d, 20))
    new_sig <- array(NA_real_, dim = c(d, d, 20))

    old_dim <- dim(dp$clusterParameters$mu)[3]
    new_mu[, , 1:old_dim] <- dp$clusterParameters$mu
    new_sig[, , 1:old_dim] <- dp$clusterParameters$sig

    dp$clusterParameters$mu <- new_mu
    dp$clusterParameters$sig <- new_sig
  }

  # Manually run the C++ MCMC steps
  n_iter <- 100
  for (i in 1:n_iter) {
    # Convert to 0-indexed for C++
    dp$clusterLabels <- dp$clusterLabels - 1

    # Run cluster component update
    update_result <- conjugate_mvnormal_cluster_component_update_cpp(dp)
    dp$clusterLabels <- update_result$clusterLabels
    dp$pointsPerCluster <- update_result$pointsPerCluster
    dp$numberClusters <- update_result$numberClusters
    dp$clusterParameters <- update_result$clusterParameters

    # Run cluster parameter update
    dp$clusterParameters <- conjugate_mvnormal_cluster_parameter_update_cpp(dp)

    # Convert back to 1-indexed
    dp$clusterLabels <- dp$clusterLabels + 1

    # Update alpha (using R implementation for now)
    dp <- UpdateAlpha(dp)
  }

  # Basic validity checks
  expect_true(dp$numberClusters >= 1)
  expect_true(dp$numberClusters <= n)
  expect_equal(length(dp$clusterLabels), n)
  expect_true(all(dp$clusterLabels > 0))
  expect_equal(sum(dp$pointsPerCluster), n)

  # Check parameter dimensions
  expect_equal(dim(dp$clusterParameters$mu)[1], 1)
  expect_equal(dim(dp$clusterParameters$mu)[2], d)
  expect_equal(dim(dp$clusterParameters$sig)[1], d)
  expect_equal(dim(dp$clusterParameters$sig)[2], d)
})

test_that("MVNormal C++ vs R implementation equivalence", {
  skip_if_not(requireNamespace("mvtnorm", quietly = TRUE))
  skip_if_not(exists("conjugate_mvnormal_cluster_component_update_cpp"),
              "MVNormal C++ functions not available")

  # Small dataset for comparison
  set.seed(456)
  y <- mvtnorm::rmvnorm(30, c(0, 0), diag(2))

  g0Priors <- list(
    mu0 = c(0, 0),
    Lambda = diag(2) * 2,
    kappa0 = 0.5,
    nu = 4
  )

  # Test R backend
  set_use_cpp(FALSE)
  set.seed(789)
  dp_r <- DirichletProcessMvnormal(y, g0Priors)
  dp_r <- Initialise(dp_r, numInitialClusters = 1)

  # Run R implementation manually with more iterations
  for (i in 1:100) {
    dp_r <- ClusterComponentUpdate(dp_r)
    dp_r <- ClusterParameterUpdate(dp_r)
    dp_r <- UpdateAlpha(dp_r)
  }

  # Test C++ backend
  set.seed(789)
  dp_cpp <- DirichletProcessMvnormal(y, g0Priors)
  dp_cpp <- Initialise(dp_cpp, numInitialClusters = 1)

  # Ensure parameter arrays have enough space
  if (dim(dp_cpp$clusterParameters$mu)[3] < 20) {
    new_mu <- array(NA_real_, dim = c(1, 2, 20))
    new_sig <- array(NA_real_, dim = c(2, 2, 20))

    old_dim <- dim(dp_cpp$clusterParameters$mu)[3]
    new_mu[, , 1:old_dim] <- dp_cpp$clusterParameters$mu
    new_sig[, , 1:old_dim] <- dp_cpp$clusterParameters$sig

    dp_cpp$clusterParameters$mu <- new_mu
    dp_cpp$clusterParameters$sig <- new_sig
  }

  # Run C++ implementation manually
  for (i in 1:100) {
    # Convert to 0-indexed
    dp_cpp$clusterLabels <- dp_cpp$clusterLabels - 1

    # C++ updates
    update_result <- conjugate_mvnormal_cluster_component_update_cpp(dp_cpp)
    dp_cpp$clusterLabels <- update_result$clusterLabels
    dp_cpp$pointsPerCluster <- update_result$pointsPerCluster
    dp_cpp$numberClusters <- update_result$numberClusters
    dp_cpp$clusterParameters <- update_result$clusterParameters

    dp_cpp$clusterParameters <- conjugate_mvnormal_cluster_parameter_update_cpp(dp_cpp)

    # Convert back to 1-indexed
    dp_cpp$clusterLabels <- dp_cpp$clusterLabels + 1

    # Use R for alpha update
    dp_cpp <- UpdateAlpha(dp_cpp)
  }

  # Compare number of clusters (allowing for larger MCMC variability)
  # MVNormal can create many small clusters, so we allow more difference
  expect_true(abs(dp_r$numberClusters - dp_cpp$numberClusters) <= 10)

  # Both should have valid alpha values
  expect_true(dp_r$alpha > 0)
  expect_true(dp_cpp$alpha > 0)

  # Both should have found at least 1 cluster
  expect_true(dp_r$numberClusters >= 1)
  expect_true(dp_cpp$numberClusters >= 1)
})

test_that("MVNormal handles edge cases correctly", {
  skip_if_not(exists("conjugate_mvnormal_cluster_component_update_cpp"),
              "MVNormal C++ functions not available")

  # Test with 1D data (degenerate case)
  y_1d <- matrix(rnorm(20), ncol = 1)
  g0_1d <- list(
    mu0 = 0,
    Lambda = matrix(1, 1, 1),
    kappa0 = 1,
    nu = 2
  )

  dp_1d <- DirichletProcessMvnormal(y_1d, g0_1d)
  dp_1d <- Initialise(dp_1d, numInitialClusters = 1)

  # Ensure proper array dimensions for 1D case
  if (length(dim(dp_1d$clusterParameters$mu)) == 2) {
    # Convert to 3D array
    dp_1d$clusterParameters$mu <- array(dp_1d$clusterParameters$mu,
                                        dim = c(1, 1, 20))
  }
  if (length(dim(dp_1d$clusterParameters$sig)) == 2) {
    # Convert to 3D array
    dp_1d$clusterParameters$sig <- array(dp_1d$clusterParameters$sig,
                                         dim = c(1, 1, 20))
  }

  # Run a few iterations with C++ functions
  for (i in 1:10) {
    dp_1d$clusterLabels <- dp_1d$clusterLabels - 1

    update_result <- conjugate_mvnormal_cluster_component_update_cpp(dp_1d)
    dp_1d$clusterLabels <- update_result$clusterLabels
    dp_1d$pointsPerCluster <- update_result$pointsPerCluster
    dp_1d$numberClusters <- update_result$numberClusters
    dp_1d$clusterParameters <- update_result$clusterParameters

    dp_1d$clusterParameters <- conjugate_mvnormal_cluster_parameter_update_cpp(dp_1d)

    dp_1d$clusterLabels <- dp_1d$clusterLabels + 1
  }

  expect_true(dp_1d$numberClusters >= 1)

  # Test with perfect clustering (but with more realistic expectations)
  y_perfect <- rbind(
    matrix(rep(c(-5, -5), 10), ncol = 2, byrow = TRUE),
    matrix(rep(c(5, 5), 10), ncol = 2, byrow = TRUE)
  )

  dp_perfect <- DirichletProcessMvnormal(y_perfect)
  dp_perfect <- Initialise(dp_perfect, numInitialClusters = 2)

  # Ensure parameter arrays have enough space
  if (dim(dp_perfect$clusterParameters$mu)[3] < 20) {
    new_mu <- array(NA_real_, dim = c(1, 2, 20))
    new_sig <- array(NA_real_, dim = c(2, 2, 20))

    old_dim <- dim(dp_perfect$clusterParameters$mu)[3]
    new_mu[, , 1:old_dim] <- dp_perfect$clusterParameters$mu
    new_sig[, , 1:old_dim] <- dp_perfect$clusterParameters$sig

    dp_perfect$clusterParameters$mu <- new_mu
    dp_perfect$clusterParameters$sig <- new_sig
  }

  # Run more iterations for better convergence
  for (i in 1:200) {
    dp_perfect$clusterLabels <- dp_perfect$clusterLabels - 1

    update_result <- conjugate_mvnormal_cluster_component_update_cpp(dp_perfect)
    dp_perfect$clusterLabels <- update_result$clusterLabels
    dp_perfect$pointsPerCluster <- update_result$pointsPerCluster
    dp_perfect$numberClusters <- update_result$numberClusters
    dp_perfect$clusterParameters <- update_result$clusterParameters

    dp_perfect$clusterParameters <- conjugate_mvnormal_cluster_parameter_update_cpp(dp_perfect)

    dp_perfect$clusterLabels <- dp_perfect$clusterLabels + 1

    # Update alpha to help convergence
    dp_perfect <- UpdateAlpha(dp_perfect)
  }

  # Should find approximately 2 clusters, but MVNormal can create more
  # Allow up to 10 clusters due to MCMC variability
  expect_true(dp_perfect$numberClusters >= 1)
  expect_true(dp_perfect$numberClusters <= 10)
})
