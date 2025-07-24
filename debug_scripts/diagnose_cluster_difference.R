# Diagnostic script to understand large cluster differences
library(dirichletprocess)
library(testthat)

cat("=== Diagnosing Large Cluster Differences ===\n")

# Use exact same data as the test
set.seed(123)
mu1 <- c(-2, -2)
mu2 <- c(2, 2)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
mvn2_size <- 25  # DEV_MODE size
test_data <- rbind(
  mvtnorm::rmvnorm(mvn2_size, mu1, sigma),
  mvtnorm::rmvnorm(mvn2_size, mu2, sigma)
)

# Single run comparison with detailed analysis
seed <- 12345
iterations <- 50  # DEV_MODE iterations

# R implementation
cat("\n=== R Implementation Analysis ===\n")
set.seed(seed)
set_use_cpp(FALSE)
dp_r <- DirichletProcessMvnormal2(test_data)
cat("R object classes:", class(dp_r), "\n")
cat("R object has m parameter:", !is.null(dp_r$m), "value:", dp_r$m, "\n")
cat("R algorithm type: Algorithm 8 (nonconjugate)\n")

# Check if R is actually using the nonconjugate path
dp_r <- Fit(dp_r, its = iterations, updatePrior = TRUE, progressBar = FALSE)
r_clusters_per_iter <- sapply(dp_r$labelsChain, function(x) length(unique(x)))
cat("R clusters per iteration:", head(r_clusters_per_iter, 10), "...\n")
cat("R final clusters:", tail(r_clusters_per_iter, 1), "\n")
cat("R mean clusters:", mean(r_clusters_per_iter), "\n")
cat("R alpha progression:", head(dp_r$alphaChain, 5), "...final:", tail(dp_r$alphaChain, 1), "\n")

# C++ implementation
cat("\n=== C++ Implementation Analysis ===\n")
set.seed(seed)
set_use_cpp(TRUE)
dp_cpp <- DirichletProcessMvnormal2(test_data)
cat("C++ can_use_cpp:", can_use_cpp(dp_cpp), "\n")
cat("C++ using_cpp:", using_cpp(), "\n")

dp_cpp <- Fit(dp_cpp, its = iterations, updatePrior = TRUE, progressBar = FALSE)
cpp_clusters_per_iter <- sapply(dp_cpp$labelsChain, function(x) length(unique(x)))
cat("C++ clusters per iteration:", head(cpp_clusters_per_iter, 10), "...\n")
cat("C++ final clusters:", tail(cpp_clusters_per_iter, 1), "\n")
cat("C++ mean clusters:", mean(cpp_clusters_per_iter), "\n")
cat("C++ alpha progression:", head(dp_cpp$alphaChain, 5), "...final:", tail(dp_cpp$alphaChain, 1), "\n")

# Analysis
cat("\n=== Difference Analysis ===\n")
r_mean <- mean(r_clusters_per_iter)
cpp_mean <- mean(cpp_clusters_per_iter)
difference <- abs(r_mean - cpp_mean)

cat("R mean clusters:", round(r_mean, 2), "\n")
cat("C++ mean clusters:", round(cpp_mean, 2), "\n")
cat("Absolute difference:", round(difference, 2), "\n")

# Check if this matches the failing test
if (difference > 10) {
  cat("❌ Large difference detected (>10) - suggests implementation issue\n")
} else {
  cat("✅ Difference within reasonable range\n")
}

# Compare convergence patterns
cat("\n=== Convergence Pattern Analysis ===\n")
cat("R cluster variance:", round(var(r_clusters_per_iter), 2), "\n")
cat("C++ cluster variance:", round(var(cpp_clusters_per_iter), 2), "\n")

# Check if either implementation is getting stuck
r_stuck <- length(unique(tail(r_clusters_per_iter, 10))) == 1
cpp_stuck <- length(unique(tail(cpp_clusters_per_iter, 10))) == 1
cat("R appears stuck in last 10 iterations:", r_stuck, "\n")
cat("C++ appears stuck in last 10 iterations:", cpp_stuck, "\n")