# Debug ClusterComponentUpdate aux parameter regeneration
library(dirichletprocess)

# Disable C++ to isolate R implementation
set_use_cpp(FALSE)
options(dirichletprocess.use_cpp_samplers = FALSE)
options(dirichletprocess.use_cpp_hierarchical = FALSE)

# Create test data
group1_data <- rbeta(5, 2, 5)
group2_data <- rbeta(5, 5, 2)
hierarchical_data <- list(group1_data, group2_data)

# Create hierarchical beta DP
dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)

cat("Before ClusterComponentUpdate:\n")
cat("Individual DP 1 aux[[1]] names:", names(dp$indDP[[1]]$aux[[1]]), "\n")
cat("aux[[1]] has mu:", "mu" %in% names(dp$indDP[[1]]$aux[[1]]), "\n")
cat("aux[[1]] has nu:", "nu" %in% names(dp$indDP[[1]]$aux[[1]]), "\n")

# Manually test PriorDraw on the mixing distribution
md <- dp$indDP[[1]]$mixingDistribution
cat("\nTesting PriorDraw on mixing distribution:\n")
test_draw <- PriorDraw(md, 1)
cat("Test draw names:", names(test_draw), "\n")

# Try ClusterComponentUpdate on individual DP
cat("\nAttempting ClusterComponentUpdate on individual DP...\n")
tryCatch({
  dp$indDP[[1]] <- ClusterComponentUpdate(dp$indDP[[1]])
  cat("Individual DP ClusterComponentUpdate succeeded\n")
  
  cat("After ClusterComponentUpdate:\n")
  cat("aux[[1]] names:", names(dp$indDP[[1]]$aux[[1]]), "\n")
  cat("aux[[1]] has mu:", "mu" %in% names(dp$indDP[[1]]$aux[[1]]), "\n")
  cat("aux[[1]] has nu:", "nu" %in% names(dp$indDP[[1]]$aux[[1]]), "\n")
  
}, error = function(e) {
  cat("Individual DP ClusterComponentUpdate failed:", e$message, "\n")
  print(e)
})

# Try ClusterComponentUpdate on hierarchical DP
cat("\nAttempting ClusterComponentUpdate on hierarchical DP...\n")
tryCatch({
  dp <- ClusterComponentUpdate(dp)
  cat("Hierarchical ClusterComponentUpdate succeeded\n")
}, error = function(e) {
  cat("Hierarchical ClusterComponentUpdate failed:", e$message, "\n")
  print(e)
})

cat("Debugging complete\n")