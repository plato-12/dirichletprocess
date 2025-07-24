# Simple consistency check for Normal distribution
library(dirichletprocess)

# Quick check with same seed
set.seed(123)
test_data <- rnorm(50)

# R version
set_use_cpp(FALSE)
options(dirichletprocess.use_cpp_samplers = FALSE)
dp_r <- DirichletProcessGaussian(test_data)
dp_r <- Fit(dp_r, its = 10, progressBar = FALSE)

# C++ version
set_use_cpp(TRUE)
enable_cpp_samplers(TRUE)
set.seed(123)
dp_cpp <- DirichletProcessGaussian(test_data)
dp_cpp <- Fit(dp_cpp, its = 10, progressBar = FALSE)

# Compare results
r_clusters <- sapply(dp_r$labelsChain, function(x) length(unique(x)))
cpp_clusters <- sapply(dp_cpp$labelsChain, function(x) length(unique(x)))

cat("R clusters:", r_clusters, "\n")
cat("C++ clusters:", cpp_clusters, "\n")
cat("Mean R clusters:", mean(r_clusters), "\n")
cat("Mean C++ clusters:", mean(cpp_clusters), "\n")
cat("Difference:", abs(mean(r_clusters) - mean(cpp_clusters)), "\n")
cat("Alpha R:", mean(dp_r$alphaChain), "\n")
cat("Alpha C++:", mean(dp_cpp$alphaChain), "\n")
cat("Alpha diff:", abs(mean(dp_r$alphaChain) - mean(dp_cpp$alphaChain)), "\n")