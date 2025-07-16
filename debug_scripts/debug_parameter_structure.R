# Debug parameter structure differences
library(dirichletprocess)

# Test 1: Simple validation case (works)
cat("=== TEST 1: Simple validation case ===\n")
test_data <- matrix(rnorm(20 * 2), ncol = 2)
prior_params <- list(
  mu0 = rep(0, 2),
  kappa0 = 1,
  nu = 3,
  Lambda = diag(2),
  covModel = "FULL"
)
md <- MvnormalCreate(prior_params)
dp <- DirichletProcessCreate(test_data, md)
dp <- Initialise(dp, numInitialClusters = 2)

cat("Original mu dimensions:", dim(dp$clusterParameters$mu), "\n")
cat("Original sig dimensions:", dim(dp$clusterParameters$sig), "\n")
cat("numberClusters:", dp$numberClusters, "\n")

# Test parameter extraction logic from LikelihoodDP
clusters_parameters <- dp$clusterParameters
cat("inherits(dp, 'mvnormal'):", inherits(dp, "mvnormal"), "\n")
cat("is.list(clusters_parameters):", is.list(clusters_parameters), "\n")

if (inherits(dp, "mvnormal") && is.list(clusters_parameters)) {
  cat("Extracting active clusters...\n")
  active_params <- list()
  for (i in seq_along(clusters_parameters)) {
    param_dims <- dim(clusters_parameters[[i]])
    cat("Parameter", i, "dimensions:", param_dims, "\n")
    if (length(param_dims) == 3 && param_dims[3] > dp$numberClusters) {
      cat("Extracting first", dp$numberClusters, "clusters\n")
      if (dp$numberClusters == 1) {
        active_params[[i]] <- array(clusters_parameters[[i]][, , 1:dp$numberClusters, drop = FALSE],
                                    dim = c(param_dims[1], param_dims[2], dp$numberClusters))
      } else {
        active_params[[i]] <- clusters_parameters[[i]][, , 1:dp$numberClusters, drop = FALSE]
      }
      cat("Extracted dimensions:", dim(active_params[[i]]), "\n")
    } else {
      active_params[[i]] <- clusters_parameters[[i]]
    }
  }
  clusters_parameters <- active_params
}

cat("Final mu dimensions:", dim(clusters_parameters$mu), "\n")
cat("Final sig dimensions:", dim(clusters_parameters$sig), "\n")

# Test likelihood calculation with extracted parameters
cat("Testing likelihood calculation with extracted parameters...\n")
set_use_cpp(FALSE)  # Use R implementation
tryCatch({
  lik <- Likelihood(dp$mixingDistribution, dp$data[1, , drop=FALSE], clusters_parameters)
  cat("Likelihood result length:", length(lik), "\n")
  cat("Likelihood values:", lik, "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\n=== TEST 2: Manual vapply test ===\n")
# Test the exact vapply call that's failing
tryCatch({
  clusters_parameters <- dp$clusterParameters
  
  # Manual extraction like in LikelihoodDP
  if (inherits(dp, "mvnormal") && is.list(clusters_parameters)) {
    active_params <- list()
    for (i in seq_along(clusters_parameters)) {
      param_dims <- dim(clusters_parameters[[i]])
      if (length(param_dims) == 3 && param_dims[3] > dp$numberClusters) {
        if (dp$numberClusters == 1) {
          active_params[[i]] <- array(clusters_parameters[[i]][, , 1:dp$numberClusters, drop = FALSE],
                                      dim = c(param_dims[1], param_dims[2], dp$numberClusters))
        } else {
          active_params[[i]] <- clusters_parameters[[i]][, , 1:dp$numberClusters, drop = FALSE]
        }
      } else {
        active_params[[i]] <- clusters_parameters[[i]]
      }
    }
    clusters_parameters <- active_params
  }
  
  # Test one call to the function used in vapply
  i <- 1
  lik <- Likelihood(dp$mixingDistribution, dp$data[i, , drop=FALSE], clusters_parameters)
  cat("Single likelihood call result length:", length(lik), "\n")
  cat("Expected length:", dp$numberClusters, "\n")
  cat("Likelihood values:", lik, "\n")
  
  # Test if the vapply call works
  likelihoodValues <- vapply(seq_len(min(3, nrow(dp$data))),
                             function(i) {
                               lik <- Likelihood(dp$mixingDistribution,
                                                 dp$data[i, , drop=FALSE],
                                                 clusters_parameters)
                               cat("Row", i, "likelihood length:", length(lik), "\n")
                               # Ensure we only return values for active clusters
                               if (length(lik) > dp$numberClusters) {
                                 lik[1:dp$numberClusters]
                               } else {
                                 lik
                               }
                             },
                             numeric(dp$numberClusters))
  cat("vapply result dimensions:", dim(likelihoodValues), "\n")
  
}, error = function(e) {
  cat("vapply Error:", e$message, "\n")
})

cat("\n=== TEST 3: LikelihoodDP calculation ===\n")
# Test 3: LikelihoodDP calculation (what actually fails)
tryCatch({
  lik_dp <- LikelihoodDP(dp)
  cat("LikelihoodDP result dimensions:", dim(lik_dp), "\n")
}, error = function(e) {
  cat("LikelihoodDP Error:", e$message, "\n")
})