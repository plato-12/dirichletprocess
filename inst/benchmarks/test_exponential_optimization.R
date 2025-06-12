# Quick test script for exponential optimizations
library(dirichletprocess)

cat("\n=== QUICK EXPONENTIAL OPTIMIZATION TEST ===\n\n")

# Test 1: Basic functionality
cat("1. Testing basic functionality...\n")
set.seed(123)
data <- rexp(100, rate = 2)

# Try C++ version
set_use_cpp(TRUE)
dp_cpp <- DirichletProcessExponential(data)
dp_cpp <- Fit(dp_cpp, 50, progressBar = FALSE)
cat("   C++ version completed successfully\n")

# Test 2: Quick performance comparison
cat("\n2. Quick performance test (n=500, 100 iterations)...\n")
set.seed(456)
data <- c(rexp(250, rate = 1), rexp(250, rate = 5))

# Time R version
set_use_cpp(FALSE)
time_r <- system.time({
  dp_r <- DirichletProcessExponential(data)
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)
})["elapsed"]

# Time C++ version
set_use_cpp(TRUE)
time_cpp <- system.time({
  dp_cpp <- DirichletProcessExponential(data)
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)
})["elapsed"]

speedup <- as.numeric(time_r) / as.numeric(time_cpp)

cat(sprintf("   R time:   %.3f seconds\n", time_r))
cat(sprintf("   C++ time: %.3f seconds\n", time_cpp))
cat(sprintf("   Speedup:  %.1fx\n", speedup))

if (speedup > 1) {
  cat("\n✓ SUCCESS: C++ is faster than R!\n")
} else {
  cat("\n✗ WARNING: C++ is still slower than R. Check compilation.\n")
}

# Test 3: Component timing
cat("\n3. Component-level timing (1000 likelihood calls)...\n")
n <- 100
data <- matrix(rexp(n, rate = 3), ncol = 1)

# Initialize both
set_use_cpp(FALSE)
dp_r <- DirichletProcessExponential(data)
dp_r <- Initialise(dp_r)

set_use_cpp(TRUE)
dp_cpp <- DirichletProcessExponential(data)
dp_cpp <- Initialise(dp_cpp)

# Time likelihood calculations
set_use_cpp(FALSE)
time_lik_r <- system.time({
  for (rep in 1:10) {
    for (i in 1:n) {
      Likelihood(dp_r$mixingDistribution,
                 data[i, , drop = FALSE],
                 dp_r$clusterParameters)
    }
  }
})["elapsed"] / 10

set_use_cpp(TRUE)
time_lik_cpp <- system.time({
  for (rep in 1:10) {
    for (i in 1:n) {
      Likelihood(dp_cpp$mixingDistribution,
                 data[i, , drop = FALSE],
                 dp_cpp$clusterParameters)
    }
  }
})["elapsed"] / 10

lik_speedup <- as.numeric(time_lik_r) / as.numeric(time_lik_cpp)
cat(sprintf("   Likelihood speedup: %.1fx\n", lik_speedup))

cat("\n=== TEST COMPLETE ===\n")
