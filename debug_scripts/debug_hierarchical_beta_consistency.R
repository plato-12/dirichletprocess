# Debug hierarchical beta R/C++ consistency issue
library(dirichletprocess)

# Create test data like in the test
group1_data <- rbeta(25, 2, 5)
group2_data <- rbeta(25, 5, 2) 
group3_data <- rbeta(25, 1, 1)
hierarchical_data <- list(group1_data, group2_data, group3_data)

cat("Testing R implementation...\n")
set_use_cpp(FALSE)
options(dirichletprocess.use_cpp_samplers = FALSE)
options(dirichletprocess.use_cpp_hierarchical = FALSE)

tryCatch({
  dp_r <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)
  cat("R DP created successfully\n")
  
  # Try one iteration
  dp_r <- Fit(dp_r, its = 1, updatePrior = TRUE)
  cat("R fitting succeeded for 1 iteration\n")
  
}, error = function(e) {
  cat("R implementation failed with error:", e$message, "\n")
  print(e)
})

cat("\nTesting C++ implementation...\n")
set_use_cpp(TRUE)
options(dirichletprocess.use_cpp_samplers = TRUE)
options(dirichletprocess.use_cpp_hierarchical = TRUE)

tryCatch({
  dp_cpp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)
  cat("C++ DP created successfully\n")
  
  # Try one iteration
  dp_cpp <- Fit(dp_cpp, its = 1, updatePrior = TRUE)
  cat("C++ fitting succeeded for 1 iteration\n")
  
}, error = function(e) {
  cat("C++ implementation failed with error:", e$message, "\n")
  print(e)
})

cat("Debugging complete\n")