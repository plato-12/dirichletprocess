# Debug with exact test conditions that are failing
library(dirichletprocess)

cat("=== Debug with Exact Test Conditions ===\n")

# Use EXACTLY the same data as the failing test
set.seed(123)
mu1 <- c(-2, -2)
mu2 <- c(2, 2)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
mvn2_size <- 25  # DEV_MODE size
test_data <- rbind(
  mvtnorm::rmvnorm(mvn2_size, mu1, sigma),
  mvtnorm::rmvnorm(mvn2_size, mu2, sigma)
)

cat("Test data dimensions:", dim(test_data), "\n")
cat("Test data sample:\n")
print(head(test_data, 3))

# Test with exact parameters from validate_r_cpp_consistency
seed <- 12346  # seed + run - 1, where run = 2
iterations <- 50  # DEV_MODE iterations

cat("\n=== C++ Implementation with Test Conditions ===\n")
set.seed(seed)
set_use_cpp(TRUE)
enable_cpp_samplers(TRUE)
dp_cpp <- DirichletProcessMvnormal2(test_data)

cat("Initial alpha:", dp_cpp$alpha, "\n")
cat("can_use_cpp:", can_use_cpp(dp_cpp), "\n")

# Fit with exact same parameters as test
dp_cpp <- Fit(dp_cpp, its = iterations, progressBar = FALSE)

# Analyze alpha behavior
cat("Alpha chain (first 10):", head(dp_cpp$alphaChain, 10), "\n")
cat("Alpha chain (last 10):", tail(dp_cpp$alphaChain, 10), "\n")
cat("Unique alpha values:", length(unique(dp_cpp$alphaChain)), "\n")
cat("Alpha range:", range(dp_cpp$alphaChain), "\n")

# Check if alpha is truly stuck
alpha_changes <- diff(dp_cpp$alphaChain) != 0
cat("Alpha changed between iterations:", sum(alpha_changes), "out of", length(alpha_changes), "\n")

if (sum(alpha_changes) == 0) {
  cat("❌ CONFIRMED: Alpha is stuck in C++ implementation\n")
  
  # Let's check the mcmc_params being passed
  cat("\n=== MCMC Parameters Debug ===\n")
  mcmc_params <- dirichletprocess:::prepare_mcmc_params(dp_cpp, iterations, TRUE)
  cat("update_concentration:", mcmc_params$update_concentration, "\n")
  cat("alpha_prior_shape:", mcmc_params$alpha_prior_shape, "\n")
  cat("alpha_prior_rate:", mcmc_params$alpha_prior_rate, "\n")
  
} else {
  cat("✅ Alpha is updating properly\n")
}

# Compare with R for same conditions
cat("\n=== R Implementation Comparison ===\n")
set.seed(seed)
set_use_cpp(FALSE)
dp_r <- DirichletProcessMvnormal2(test_data)
dp_r <- Fit(dp_r, its = iterations, progressBar = FALSE)
cat("R alpha range:", range(dp_r$alphaChain), "\n")
cat("R unique alpha values:", length(unique(dp_r$alphaChain)), "\n")