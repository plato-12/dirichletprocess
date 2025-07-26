# MCMC Chain Storage Diagnostic
# Investigate why chains are not being stored in basic tests

library(dirichletprocess)

# Test 1: Simple Gaussian DP
cat("=== Testing Basic Gaussian DP ===\n")
data_test <- rnorm(10)
dpobj <- DirichletProcessGaussian(data_test)

cat("Before Fit:\n")
cat("  alphaChain length:", length(dpobj$alphaChain), "\n")
cat("  weightsChain length:", length(dpobj$weightsChain), "\n")
cat("  clusterParametersChain length:", length(dpobj$clusterParametersChain), "\n")

# Fit with 10 iterations
dpobj <- Fit(dpobj, its = 10, updatePrior = FALSE, verbose = TRUE)

cat("After Fit (10 iterations):\n")
cat("  alphaChain length:", length(dpobj$alphaChain), "\n")
cat("  weightsChain length:", length(dpobj$weightsChain), "\n")
cat("  clusterParametersChain length:", length(dpobj$clusterParametersChain), "\n")
cat("  Class:", class(dpobj), "\n")

# Check if chains contain actual data
if (length(dpobj$alphaChain) > 0) {
  cat("  alphaChain first 5 values:", head(dpobj$alphaChain, 5), "\n")
}

if (length(dpobj$clusterParametersChain) > 0) {
  cat("  clusterParametersChain structure:\n")
  str(dpobj$clusterParametersChain[1:min(3, length(dpobj$clusterParametersChain))])
}

# Test 2: Check C++ vs R behavior
cat("\n=== Testing C++ vs R Implementation ===\n")

# Test R implementation
cat("--- R Implementation ---\n")
set_use_cpp(FALSE)
dpobj_r <- DirichletProcessGaussian(data_test)
dpobj_r <- Fit(dpobj_r, its = 5, updatePrior = FALSE, verbose = FALSE)
cat("R implementation chains:\n")
cat("  alphaChain length:", length(dpobj_r$alphaChain), "\n")
cat("  clusterParametersChain length:", length(dpobj_r$clusterParametersChain), "\n")

# Test C++ implementation
cat("--- C++ Implementation ---\n")
set_use_cpp(TRUE)
dpobj_cpp <- DirichletProcessGaussian(data_test)
dpobj_cpp <- Fit(dpobj_cpp, its = 5, updatePrior = FALSE, verbose = FALSE)
cat("C++ implementation chains:\n")
cat("  alphaChain length:", length(dpobj_cpp$alphaChain), "\n")
cat("  clusterParametersChain length:", length(dpobj_cpp$clusterParametersChain), "\n")

# Test 3: Check what using_cpp() returns
cat("\n=== C++ Status ===\n")
cat("using_cpp():", using_cpp(), "\n")
cat("C++ available:", can_use_cpp(), "\n")

# Test 4: Weibull (non-conjugate) behavior
cat("\n=== Testing Weibull (Non-conjugate) ===\n")
data_weibull <- rweibull(10, 1, 1)
priorParameters <- matrix(c(1,1,1), ncol=3)

# R implementation
set_use_cpp(FALSE)
dpobj_weibull_r <- DirichletProcessWeibull(data_weibull, g0Priors = c(1,1,1))
dpobj_weibull_r <- Fit(dpobj_weibull_r, its = 5, updatePrior = TRUE, verbose = FALSE)
cat("Weibull R implementation:\n")
cat("  alphaChain length:", length(dpobj_weibull_r$alphaChain), "\n")
cat("  clusterParametersChain length:", length(dpobj_weibull_r$clusterParametersChain), "\n")
if (!is.null(dpobj_weibull_r$priorParametersChain)) {
  cat("  priorParametersChain length:", length(dpobj_weibull_r$priorParametersChain), "\n")
}

# C++ implementation  
set_use_cpp(TRUE)
dpobj_weibull_cpp <- DirichletProcessWeibull(data_weibull, g0Priors = c(1,1,1))
dpobj_weibull_cpp <- Fit(dpobj_weibull_cpp, its = 5, updatePrior = TRUE, verbose = FALSE)
cat("Weibull C++ implementation:\n")
cat("  alphaChain length:", length(dpobj_weibull_cpp$alphaChain), "\n")
cat("  clusterParametersChain length:", length(dpobj_weibull_cpp$clusterParametersChain), "\n")
if (!is.null(dpobj_weibull_cpp$priorParametersChain)) {
  cat("  priorParametersChain length:", length(dpobj_weibull_cpp$priorParametersChain), "\n")
}

cat("\n=== Diagnostic Complete ===\n")