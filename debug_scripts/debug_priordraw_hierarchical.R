# Debug PriorDraw.hierarchical structure
library(dirichletprocess)

# Create simple test data
group1_data <- rbeta(5, 2, 5)
group2_data <- rbeta(5, 5, 2)
hierarchical_data <- list(group1_data, group2_data)

# Create hierarchical beta DP
dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)

# Check the mixing distribution structure
md <- dp$indDP[[1]]$mixingDistribution
cat("Mixing distribution class:", class(md), "\n")
cat("theta_k names:", names(md$theta_k), "\n")
cat("theta_k structure:\n")
str(md$theta_k)

# Test PriorDraw.hierarchical
cat("\nTesting PriorDraw.hierarchical...\n")
prior_draw <- PriorDraw(md, 1)
cat("PriorDraw result names:", names(prior_draw), "\n")
cat("PriorDraw result structure:\n")
str(prior_draw)

# Check if it has mu and nu
if (is.list(prior_draw)) {
  cat("Has 'mu':", "mu" %in% names(prior_draw), "\n")
  cat("Has 'nu':", "nu" %in% names(prior_draw), "\n")
}

# Test Likelihood.beta with this structure
if ("mu" %in% names(prior_draw) && "nu" %in% names(prior_draw)) {
  cat("Testing Likelihood.beta...\n")
  tryCatch({
    lik <- Likelihood(md, matrix(0.5, nrow=1), prior_draw)
    cat("Likelihood succeeded:", lik, "\n")
  }, error = function(e) {
    cat("Likelihood failed:", e$message, "\n")
  })
} else {
  cat("Cannot test Likelihood - missing mu/nu\n")
}