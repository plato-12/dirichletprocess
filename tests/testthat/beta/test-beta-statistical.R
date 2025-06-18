context("Beta DP Statistical Properties")

test_that("Beta DP posterior is consistent", {
  set.seed(5678)
  # Generate data with known properties
  true_a <- 4
  true_b <- 6
  n_data <- 200
  data <- rbeta(n_data, true_a, true_b)

  # Fit model
  dp <- DirichletProcessBeta(data, verbose = FALSE)
  dp <- Fit(dp, its = 1000, progressBar = FALSE)

  # Extract posterior samples
  n_samples <- length(dp$alphaChain)
  mu_posterior <- numeric(n_samples)
  tau_posterior <- numeric(n_samples)

  # Get weighted average of parameters
  for (i in 1:n_samples) {
    if (i <= length(dp$labelsChain)) {
      labels <- dp$labelsChain[[i]]
      n_clusters <- length(unique(labels))
      weights <- table(labels) / length(labels)

      # Weighted average of mu
      mu_clusters <- dp$clusterParameters[[1]][1:n_clusters]
      mu_posterior[i] <- sum(weights * mu_clusters)

      # Weighted average of tau
      tau_clusters <- dp$clusterParameters[[2]][1:n_clusters]
      tau_posterior[i] <- sum(weights * tau_clusters)
    }
  }

  # Check posterior mean
  true_mean <- true_a / (true_a + true_b)
  posterior_mean <- mean(mu_posterior[-(1:100)])  # Remove burn-in
  expect_equal(posterior_mean, true_mean, tolerance = 0.05)
})

test_that("Beta DP uncertainty quantification is correct", {
  set.seed(6789)

  # Small sample - high uncertainty
  small_data <- rbeta(10, 3, 3)
  dp_small <- DirichletProcessBeta(small_data, verbose = FALSE)
  dp_small <- Fit(dp_small, its = 500, progressBar = FALSE)

  # Large sample - low uncertainty
  large_data <- rbeta(500, 3, 3)
  dp_large <- DirichletProcessBeta(large_data, verbose = FALSE)
  dp_large <- Fit(dp_large, its = 500, progressBar = FALSE)

  # Extract parameter chains
  get_param_variance <- function(dp) {
    n_iter <- length(dp$labelsChain)
    mu_values <- numeric(n_iter)

    for (i in 1:n_iter) {
      if (dp$numberClusters == 1) {
        mu_values[i] <- dp$clusterParameters[[1]][1]
      } else {
        # Weighted average
        weights <- dp$pointsPerCluster / dp$n
        mu_values[i] <- sum(weights * dp$clusterParameters[[1]][1:dp$numberClusters])
      }
    }

    var(mu_values[-(1:100)])  # Remove burn-in
  }

  var_small <- get_param_variance(dp_small)
  var_large <- get_param_variance(dp_large)

  # Variance should decrease with more data
  expect_true(var_large < var_small / 2)
})

test_that("Beta DP prior sensitivity analysis", {
  set.seed(7890)
  data <- rbeta(50, 4, 6)

  # Different priors
  priors <- list(
    weak = c(1, 1),      # Weak prior
    moderate = c(2, 8),  # Moderate prior
    strong = c(10, 40)   # Strong prior
  )

  results <- list()

  for (prior_name in names(priors)) {
    dp <- DirichletProcessBeta(
      data,
      g0Priors = priors[[prior_name]],
      verbose = FALSE
    )
    dp <- Fit(dp, its = 500, progressBar = FALSE)

    # Get posterior mean
    if (dp$numberClusters == 1) {
      post_mean <- dp$clusterParameters[[1]][1]
    } else {
      weights <- dp$pointsPerCluster / dp$n
      post_mean <- sum(weights * dp$clusterParameters[[1]][1:dp$numberClusters])
    }

    results[[prior_name]] <- post_mean
  }

  # Strong prior should pull estimate toward prior mean
  prior_mean_strong <- priors$strong[1] / sum(priors$strong)
  expect_true(abs(results$strong - prior_mean_strong) <
                abs(results$weak - prior_mean_strong))
})
