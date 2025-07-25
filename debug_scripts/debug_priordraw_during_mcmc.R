# Debug PriorDraw behavior during MCMC
library(dirichletprocess)
source("tests/testthat/helper-testing.R")

# Disable C++ to isolate R implementation
set_use_cpp(FALSE)
options(dirichletprocess.use_cpp_samplers = FALSE)
options(dirichletprocess.use_cpp_hierarchical = FALSE)

# Create test data exactly like in the test
set.seed(123)
group1_data <- rbeta(25, 2, 5)
group2_data <- rbeta(25, 5, 2) 
group3_data <- rbeta(25, 1, 1)
hierarchical_data <- list(group1_data, group2_data, group3_data)
attr(hierarchical_data, "maxY") <- 1

dp_r <- create_dp_object("hierarchical_beta", hierarchical_data)

# Run first iteration successfully
dp_r <- Fit(dp_r, its = 1, updatePrior = TRUE)
cat("After iteration 1 - aux still has names\n")

# Test PriorDraw directly on the mixing distributions after iteration 1
for (i in 1:3) {
  cat("Testing PriorDraw on DP", i, "mixing distribution:\n")
  md <- dp_r$indDP[[i]]$mixingDistribution
  cat("  theta_k names:", names(md$theta_k), "\n")
  
  test_draw <- PriorDraw(md, 1)
  cat("  PriorDraw result names:", names(test_draw), "\n")
  
  # Check if theta_k structure changed
  if (is.null(names(md$theta_k))) {
    cat("  WARNING: theta_k lost its names!\n")
  }
  
  cat("\n")
}

# Run second iteration (this is where it breaks)
cat("Starting iteration 2...\n")
tryCatch({
  dp_r <- Fit(dp_r, its = 1, updatePrior = TRUE)
  cat("Iteration 2 succeeded\n")
}, error = function(e) {
  cat("Iteration 2 failed:", e$message, "\n")
})

# Check PriorDraw after iteration 2
for (i in 1:3) {
  cat("After iteration 2 - Testing PriorDraw on DP", i, ":\n")
  md <- dp_r$indDP[[i]]$mixingDistribution
  cat("  theta_k names:", names(md$theta_k), "\n")
  
  test_draw <- PriorDraw(md, 1)
  cat("  PriorDraw result names:", names(test_draw), "\n")
  cat("\n")
}

cat("Debugging complete\n")