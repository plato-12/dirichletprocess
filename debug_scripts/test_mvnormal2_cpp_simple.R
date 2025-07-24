# Simple MVNormal2 C++ functionality test
library(dirichletprocess)

# Test data
set.seed(123)
mu1 <- c(-1, -1)
mu2 <- c(1, 1)
sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
test_data <- rbind(
  mvtnorm::rmvnorm(10, mu1, sigma),
  mvtnorm::rmvnorm(10, mu2, sigma)
)

cat("=== Testing MVNormal2 R vs C++ Basic Functionality ===\n")

# Test R implementation
cat("\n1. Testing R implementation...\n")
set_use_cpp(FALSE)
dp_r <- DirichletProcessMvnormal2(test_data)
cat("R DP object created successfully\n")
dp_r <- Fit(dp_r, its = 5, progressBar = FALSE)
cat("R fit completed, clusters:", length(unique(dp_r$clusterLabels)), "\n")

# Test C++ implementation
cat("\n2. Testing C++ implementation...\n")
set_use_cpp(TRUE)
dp_cpp <- DirichletProcessMvnormal2(test_data)
cat("C++ DP object created successfully\n")
cat("can_use_cpp:", can_use_cpp(dp_cpp), "\n")
cat("using_cpp:", using_cpp(), "\n")

dp_cpp <- Fit(dp_cpp, its = 5, progressBar = FALSE)
cat("C++ fit completed, clusters:", length(unique(dp_cpp$clusterLabels)), "\n")

# Compare basic results
cat("\n3. Basic comparison:\n")
cat("R final alpha:", tail(dp_r$alphaChain, 1), "\n")
cat("C++ final alpha:", tail(dp_cpp$alphaChain, 1), "\n")
cat("R final clusters:", length(unique(dp_r$clusterLabels)), "\n")
cat("C++ final clusters:", length(unique(dp_cpp$clusterLabels)), "\n")

cat("\n=== Test completed successfully ===\n")