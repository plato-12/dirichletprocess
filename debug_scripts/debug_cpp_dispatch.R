# Debug C++ dispatch for mvnormal
library(dirichletprocess)

# Enable C++
set_use_cpp(TRUE)
enable_cpp_samplers()

# Create test data
set.seed(42)
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

cat("=== TESTING DIRECT C++ DISPATCH ===\n")
cat("Class of mdObj:", class(md), "\n")
cat("mu dimensions:", dim(dp$clusterParameters$mu), "\n")
cat("sig dimensions:", dim(dp$clusterParameters$sig), "\n")

# Extract parameters like LikelihoodDP does
clusters_parameters <- dp$clusterParameters
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
  names(clusters_parameters) <- names(dp$clusterParameters)
}

cat("After extraction:\n")
cat("mu dimensions:", dim(clusters_parameters$mu), "\n")
cat("sig dimensions:", dim(clusters_parameters$sig), "\n")

# Test direct call to Likelihood function
cat("Testing direct Likelihood call...\n")
tryCatch({
  result <- Likelihood(md, dp$data[1, , drop=FALSE], clusters_parameters)
  cat("Direct call result length:", length(result), "\n")
  cat("Direct call result:", result, "\n")
}, error = function(e) {
  cat("Direct call error:", e$message, "\n")
  cat("Call stack:\n")
  print(sys.calls())
})

# Test the C++ dispatch manually
cat("\n=== TESTING MANUAL C++ DISPATCH ===\n")
x <- as.numeric(dp$data[1, ])
dist_type <- class(md)[class(md) != "list" & class(md) != "MixingDistribution"][1]
cat("dist_type:", dist_type, "\n")
cat("x:", x, "\n")
cat("theta structure:\n")
str(clusters_parameters)

# Check if the mvnormal handler is being called
if (dist_type == "mvnormal" || any(grepl("mvnormal", class(md)))) {
  cat("mvnormal handler would be called\n")
  
  # Test the parameter extraction
  mu_array <- clusters_parameters$mu
  sig_array <- clusters_parameters$sig
  mu_dim <- dim(mu_array)
  
  cat("mu_dim:", mu_dim, "\n")
  cat("sig_dim:", dim(sig_array), "\n")
  
  # Extract number of clusters
  if (is.null(mu_dim) || length(mu_dim) < 3) {
    if (is.null(mu_dim)) {
      num_clusters <- 1
    } else if (length(mu_dim) == 2) {
      num_clusters <- mu_dim[2]
    } else {
      num_clusters <- 1
    }
  } else {
    num_clusters <- mu_dim[3]
  }
  
  cat("num_clusters:", num_clusters, "\n")
  
  # Test getting the C++ function
  tryCatch({
    pkg_ns <- getNamespace("dirichletprocess")
    mvnormal_likelihood_cpp <- get("mvnormal_likelihood_cpp", pkg_ns)
    cat("C++ function obtained successfully\n")
    
    # Test calling it for one cluster
    k <- 1
    if (length(mu_dim) == 3) {
      cluster_mu <- mu_array[, , k]
      cluster_sig <- sig_array[, , k]
    } else {
      cluster_mu <- mu_array[, k]
      cluster_sig <- sig_array[, , k]
    }
    
    cat("cluster_mu:", cluster_mu, "\n")
    cat("cluster_sig dimensions:", dim(cluster_sig), "\n")
    
    # Test C++ call
    x_mat <- matrix(x, nrow = 1)
    mu_vec <- as.vector(cluster_mu)
    sig_mat <- cluster_sig
    
    cat("Testing C++ call...\n")
    result <- mvnormal_likelihood_cpp(x_mat, mu_vec, sig_mat)
    cat("C++ call result:", result, "\n")
    
  }, error = function(e) {
    cat("C++ call error:", e$message, "\n")
  })
}