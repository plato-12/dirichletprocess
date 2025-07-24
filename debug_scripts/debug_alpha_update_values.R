# Debug the exact values used in C++ alpha update
library(dirichletprocess)

cat("=== Debug Alpha Update Numerical Values ===\n")

# Recreate the exact problematic scenario
set.seed(123)
mu1 <- c(-2, -2)
mu2 <- c(2, 2)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
mvn2_size <- 25
test_data <- rbind(
  mvtnorm::rmvnorm(mvn2_size, mu1, sigma),
  mvtnorm::rmvnorm(mvn2_size, mu2, sigma)
)

# Set up C++ run
seed <- 12346
set.seed(seed)
set_use_cpp(TRUE)
dp_cpp <- DirichletProcessMvnormal2(test_data)

cat("Data dimensions:", dim(test_data), "\n")
cat("Initial alpha:", dp_cpp$alpha, "\n")
cat("Alpha prior shape:", dp_cpp$alphaPriorParameters[1], "\n")
cat("Alpha prior rate:", dp_cpp$alphaPriorParameters[2], "\n")

# Run one iteration to see what's happening
dp_cpp_1iter <- Fit(dp_cpp, its = 1, progressBar = FALSE)
cat("Alpha after 1 iteration:", tail(dp_cpp_1iter$alphaChain, 1), "\n")

# Let's manually calculate what the alpha update SHOULD produce
# Using the same logic as the C++ code
current_alpha <- dp_cpp$alpha
n_rows <- nrow(test_data)
alpha_prior_shape <- dp_cpp$alphaPriorParameters[1]  # 2
alpha_prior_rate <- dp_cpp$alphaPriorParameters[2]   # 4

cat("\n=== Manual Alpha Update Calculation ===\n")
cat("Current alpha:", current_alpha, "\n")
cat("Data rows:", n_rows, "\n")
cat("Prior shape:", alpha_prior_shape, "\n")
cat("Prior rate:", alpha_prior_rate, "\n")

# Estimate number of clusters (C++ might be using a different count)
# Let's see how many clusters the C++ thinks it has after 1 iteration
clusters_after_1 <- length(unique(dp_cpp_1iter$clusterLabels))
cat("Clusters after 1 iteration:", clusters_after_1, "\n")

# Test the beta sampling step
set.seed(seed + 100)  # Different seed for this test
x_sample <- rbeta(1, current_alpha + 1.0, n_rows)
cat("Beta sample x:", x_sample, "\n")
log_x <- log(x_sample)
cat("log(x):", log_x, "\n")

# Calculate the probabilities
pi1 <- alpha_prior_shape + clusters_after_1 - 1.0
pi2 <- n_rows * (alpha_prior_rate - log_x)
pi_ratio <- pi1 / (pi1 + pi2)

cat("pi1:", pi1, "\n")
cat("pi2:", pi2, "\n")
cat("pi_ratio:", pi_ratio, "\n")

# Check for numerical issues
if (!is.finite(pi_ratio)) {
  cat("❌ pi_ratio is not finite!\n")
} else if (pi_ratio < 0 || pi_ratio > 1) {
  cat("❌ pi_ratio is out of bounds:", pi_ratio, "\n")
} else {
  cat("✅ pi_ratio looks OK\n")
}

# Test the gamma sampling
post_rate <- alpha_prior_rate - log_x
if (post_rate <= 0.0) {
  cat("❌ post_rate is non-positive:", post_rate, "\n")
  post_rate <- 0.001
  cat("   Corrected to:", post_rate, "\n")
}

# Sample new alpha
set.seed(seed + 200)
uniform_sample <- runif(1)
if (uniform_sample < pi_ratio) {
  post_shape <- alpha_prior_shape + clusters_after_1
} else {
  post_shape <- alpha_prior_shape + clusters_after_1 - 1.0
}

cat("post_shape:", post_shape, "\n")
cat("post_rate:", post_rate, "\n")

new_alpha <- rgamma(1, post_shape, rate = post_rate)
cat("Manually calculated new alpha:", new_alpha, "\n")

if (abs(new_alpha - current_alpha) < 1e-10) {
  cat("❌ New alpha is essentially identical to old alpha\n")
} else {
  cat("✅ New alpha is different from old alpha\n")
}