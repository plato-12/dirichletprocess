# Debug the C++ alpha update issue
library(dirichletprocess)

cat("=== C++ Alpha Update Debug ===\n")

# Simple test case
set.seed(123)
test_data <- rbind(c(1,1), c(2,2), c(3,3))

set_use_cpp(TRUE)
dp <- DirichletProcessMvnormal2(test_data)

cat("Initial alpha:", dp$alpha, "\n")
cat("Alpha prior parameters:", dp$alphaPriorParameters, "\n")
cat("updatePrior setting: TRUE (should enable alpha updates)\n")

# Test with very short run to see if alpha updates at all
dp_fitted <- Fit(dp, its = 5, updatePrior = TRUE, progressBar = FALSE)

cat("\n=== Alpha Chain Analysis ===\n")
print(dp_fitted$alphaChain)
cat("Alpha changed from initial:", dp$alpha != dp_fitted$alphaChain[1], "\n")
cat("Alpha varies across chain:", length(unique(dp_fitted$alphaChain)) > 1, "\n")

# Check if the issue is with the updatePrior parameter
cat("\n=== Testing with updatePrior = FALSE ===\n")
set.seed(123)
dp2 <- DirichletProcessMvnormal2(test_data)
dp2_fitted <- Fit(dp2, its = 5, updatePrior = FALSE, progressBar = FALSE)
cat("Alpha chain with updatePrior=FALSE:", dp2_fitted$alphaChain, "\n")
cat("Alpha varies with updatePrior=FALSE:", length(unique(dp2_fitted$alphaChain)) > 1, "\n")

# Check the R implementation for comparison
cat("\n=== R Implementation Comparison ===\n")
set_use_cpp(FALSE)
set.seed(123)
dp_r <- DirichletProcessMvnormal2(test_data)
dp_r_fitted <- Fit(dp_r, its = 5, updatePrior = TRUE, progressBar = FALSE)
cat("R alpha chain:", dp_r_fitted$alphaChain, "\n")
cat("R alpha varies:", length(unique(dp_r_fitted$alphaChain)) > 1, "\n")