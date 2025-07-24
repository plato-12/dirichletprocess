# Debug script for MVNormal NaN correlation issue
# debug_scripts/debug_mvnormal_nan_correlation.R

library(dirichletprocess)

# Reproduce the exact test setup
set.seed(123)
mu1 <- c(0, 0)
mu2 <- c(3, 3)
sigma <- diag(2)
DEV_MODE <- TRUE
mvn_size <- if (DEV_MODE) 25 else 50
BASE_ITERATIONS <- if (DEV_MODE) 50 else 200

test_data <- rbind(
  mvtnorm::rmvnorm(mvn_size, mu1, sigma),
  mvtnorm::rmvnorm(mvn_size, mu2, sigma)
)

cat("Test data dimensions:", dim(test_data), "\n")
cat("BASE_ITERATIONS:", BASE_ITERATIONS, "\n")

# Test R implementation
cat("=== Testing R Implementation ===\n")
set.seed(123)
set_use_cpp(FALSE)
options(dirichletprocess.use_cpp_samplers = FALSE)
options(dirichletprocess.use_cpp_hierarchical = FALSE)

dp_r <- DirichletProcessMvnormal(test_data)
cat("R DP object created successfully\n")
cat("Initial alpha:", dp_r$alpha, "\n")

dp_r <- Fit(dp_r, its = BASE_ITERATIONS)
cat("R MCMC completed\n")
cat("R likelihood chain length:", length(dp_r$likelihoodChain), "\n")
cat("R likelihood chain sample (first 5):", head(dp_r$likelihoodChain, 5), "\n")
cat("R likelihood chain sample (last 5):", tail(dp_r$likelihoodChain, 5), "\n")
cat("R likelihood chain has NA values:", any(is.na(dp_r$likelihoodChain)), "\n")
cat("R likelihood chain variance:", var(dp_r$likelihoodChain, na.rm = TRUE), "\n")

# Test C++ implementation
cat("\n=== Testing C++ Implementation ===\n")
set.seed(123)
set_use_cpp(TRUE)
enable_cpp_samplers(TRUE)
enable_cpp_hierarchical_samplers(TRUE)

dp_cpp <- DirichletProcessMvnormal(test_data)
cat("C++ DP object created successfully\n")
cat("Initial alpha:", dp_cpp$alpha, "\n")
cat("Using C++:", using_cpp(), "\n")

dp_cpp <- Fit(dp_cpp, its = BASE_ITERATIONS)
cat("C++ MCMC completed\n")
cat("C++ likelihood chain length:", length(dp_cpp$likelihoodChain), "\n") 
cat("C++ likelihood chain sample (first 5):", head(dp_cpp$likelihoodChain, 5), "\n")
cat("C++ likelihood chain sample (last 5):", tail(dp_cpp$likelihoodChain, 5), "\n")
cat("C++ likelihood chain has NA values:", any(is.na(dp_cpp$likelihoodChain)), "\n")
cat("C++ likelihood chain variance:", var(dp_cpp$likelihoodChain, na.rm = TRUE), "\n")

# Test correlation calculation
cat("\n=== Correlation Analysis ===\n")
cat("Both chains same length:", length(dp_r$likelihoodChain) == length(dp_cpp$likelihoodChain), "\n")

if (length(dp_r$likelihoodChain) > 0 && length(dp_cpp$likelihoodChain) > 0) {
  cat("Attempting correlation calculation...\n")
  correlation_result <- tryCatch({
    cor_val <- cor(dp_r$likelihoodChain, dp_cpp$likelihoodChain)
    cat("Correlation result:", cor_val, "\n")
    cor_val
  }, error = function(e) {
    cat("Error in correlation:", e$message, "\n")
    return(NA)
  })
  
  # Additional diagnostics
  cat("R chain all identical:", length(unique(dp_r$likelihoodChain)) == 1, "\n")
  cat("C++ chain all identical:", length(unique(dp_cpp$likelihoodChain)) == 1, "\n")
  
  if (length(unique(dp_r$likelihoodChain)) == 1) {
    cat("R chain unique value:", unique(dp_r$likelihoodChain), "\n")
  }
  if (length(unique(dp_cpp$likelihoodChain)) == 1) {
    cat("C++ chain unique value:", unique(dp_cpp$likelihoodChain), "\n")
  }
} else {
  cat("One or both chains are empty!\n")
}

cat("\n=== Additional Diagnostics ===\n")
cat("R alpha chain length:", length(dp_r$alphaChain), "\n")
cat("C++ alpha chain length:", length(dp_cpp$alphaChain), "\n")
cat("R final alpha:", tail(dp_r$alphaChain, 1), "\n")
cat("C++ final alpha:", tail(dp_cpp$alphaChain, 1), "\n")