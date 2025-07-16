# Debug Initialise step by step
library(dirichletprocess)
set_use_cpp(TRUE)

cat("=== DEBUGGING INITIALISE STEP BY STEP ===\n")

# Create E model
params_1d <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1, 1, 1),
  covModel = "E"
)

md <- MvnormalCreate(params_1d)
data_1d <- matrix(rnorm(10), ncol = 1)
dp <- DirichletProcessCreate(data_1d, md)

cat("Created DP object, calling Initialise manually...\n")

# Manual Initialise.conjugate step by step
numInitialClusters <- 2
dp$clusterLabels <- rep_len(seq_len(numInitialClusters), length.out = dp$n)
dp$numberClusters <- numInitialClusters
dp$pointsPerCluster <- vapply(seq_len(numInitialClusters), function(x) sum(dp$clusterLabels == x), numeric(1))

cat("Set cluster labels\n")

# Get parameters from prior
dp$clusterParameters <- PriorDraw(dp$mixingDistribution, numInitialClusters)

cat("Got prior parameters:\n")
print(dp$clusterParameters)

# Now the dimension checking part
mu_dim <- dim(dp$clusterParameters$mu)
sig_dim <- dim(dp$clusterParameters$sig)

cat(sprintf("mu_dim: %s\n", paste(mu_dim, collapse = "x")))
cat(sprintf("sig_dim: %s\n", paste(sig_dim, collapse = "x")))

if (is.null(mu_dim) || length(mu_dim) < 3) {
  cat("Handling dimension conversion...\n")
  
  if (is.null(mu_dim)) {
    # It's a vector, convert to proper 3D array
    d <- 1
    n_clusters <- length(dp$clusterParameters$mu)
    
    cat(sprintf("Converting vector mu to 3D array: d=%d, n_clusters=%d\n", d, n_clusters))
    
    dp$clusterParameters$mu <- array(dp$clusterParameters$mu, dim = c(1, d, n_clusters))
    
    cat("Converted mu, sig already looks correct:\n")
    print(dim(dp$clusterParameters$sig))
  }
  
  # Update dimensions
  mu_dim <- dim(dp$clusterParameters$mu)
  sig_dim <- dim(dp$clusterParameters$sig)
  
  cat(sprintf("Updated mu_dim: %s\n", paste(mu_dim, collapse = "x")))
  cat(sprintf("Updated sig_dim: %s\n", paste(sig_dim, collapse = "x")))
}

# Now test the expansion logic
min_slots <- max(50, dp$n, numInitialClusters * 10)
current_clusters <- mu_dim[3]

cat(sprintf("min_slots: %d, current_clusters: %d\n", min_slots, current_clusters))

if (current_clusters < min_slots) {
  cat("Need to expand arrays...\n")
  
  # Handle dimension access safely
  if (is.null(mu_dim) || length(mu_dim) < 2) {
    d <- ncol(dp$data)
  } else {
    d <- mu_dim[2]
  }
  
  cat(sprintf("d = %d\n", d))
  
  # Create new arrays with more space
  new_mu <- array(NA_real_, dim = c(1, d, min_slots))
  
  # For E model, this should be constrained
  if (exists("priorParameters", dp$mixingDistribution) && 
      !is.null(dp$mixingDistribution$priorParameters$covModel) &&
      dp$mixingDistribution$priorParameters$covModel != "FULL") {
    
    cat("E model - using constrained path\n")
    
    # Get number of parameters for this covariance model
    nParams <- dim(dp$clusterParameters$sig)[1]
    new_sig <- array(NA_real_, dim = c(nParams, min_slots))
    
    cat(sprintf("nParams: %d, creating sig array: %dx%d\n", nParams, nParams, min_slots))
    
    # Copy existing parameters
    cat("Copying existing mu...\n")
    new_mu[, , 1:current_clusters] <- dp$clusterParameters$mu
    
    cat("Copying existing sig...\n")
    print(dim(dp$clusterParameters$sig))
    print(dim(new_sig))
    
    # This is the problematic line
    cat(sprintf("Trying to copy sig - sig_dim: %s\n", paste(sig_dim, collapse = "x")))
    
    if (length(sig_dim) == 3) {
      # Convert 3D sig to 2D for constrained models
      sig_2d <- matrix(dp$clusterParameters$sig, nrow = sig_dim[1], ncol = sig_dim[3])
      new_sig[, 1:sig_dim[3]] <- sig_2d
    } else {
      new_sig[, 1:sig_dim[2]] <- dp$clusterParameters$sig
    }
    
    cat("Expansion complete\n")
  }
}

cat("\n=== STEP BY STEP DEBUG COMPLETE ===\n")