# Test Algorithm Selection in Modified C++ Implementation
library(dirichletprocess)

cat("=== Testing Algorithm Selection in C++ Implementation ===\n")

# Test 1: Normal Distribution (Conjugate - should use Algorithm 4)
cat("\n--- Test 1: Normal Distribution (Conjugate) ---\n")
set.seed(123)
test_data <- rnorm(50)

set_use_cpp(TRUE)
enable_cpp_samplers(TRUE)

dp_normal <- DirichletProcessGaussian(test_data)
cat("Distribution type:", class(dp_normal)[1], "\n")
cat("Conjugate:", "conjugate" %in% class(dp_normal), "\n")

# Run a short MCMC chain
dp_normal_fitted <- Fit(dp_normal, its = 20, progressBar = FALSE)
normal_clusters <- sapply(dp_normal_fitted$labelsChain, function(x) length(unique(x)))
cat("Cluster counts:", normal_clusters, "\n")
cat("Mean clusters:", mean(normal_clusters), "\n")

# Test 2: Beta Distribution (Non-conjugate - should use Algorithm 8)
cat("\n--- Test 2: Beta Distribution (Non-conjugate) ---\n")
set.seed(123)
test_data_beta <- rbeta(50, 2, 3)

dp_beta <- DirichletProcessBeta(test_data_beta)
cat("Distribution type:", class(dp_beta)[1], "\n")
cat("Conjugate:", "conjugate" %in% class(dp_beta), "\n")

# Run a short MCMC chain
dp_beta_fitted <- Fit(dp_beta, its = 20, progressBar = FALSE)
beta_clusters <- sapply(dp_beta_fitted$labelsChain, function(x) length(unique(x)))
cat("Cluster counts:", beta_clusters, "\n")
cat("Mean clusters:", mean(beta_clusters), "\n")

# Test 3: Exponential Distribution (Conjugate - should use Algorithm 4)
cat("\n--- Test 3: Exponential Distribution (Conjugate) ---\n")
set.seed(123)
test_data_exp <- rexp(50, rate = 1)

dp_exp <- DirichletProcessExponential(test_data_exp)
cat("Distribution type:", class(dp_exp)[1], "\n")
cat("Conjugate:", "conjugate" %in% class(dp_exp), "\n")

# Run a short MCMC chain
dp_exp_fitted <- Fit(dp_exp, its = 20, progressBar = FALSE)
exp_clusters <- sapply(dp_exp_fitted$labelsChain, function(x) length(unique(x)))
cat("Cluster counts:", exp_clusters, "\n")
cat("Mean clusters:", mean(exp_clusters), "\n")

cat("\n=== Algorithm Selection Test Complete ===\n")
cat("If successful, conjugate distributions (Normal, Exponential) now use Algorithm 4\n")
cat("and non-conjugate distributions (Beta) use Algorithm 8.\n")