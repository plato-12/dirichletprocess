# Debug create_dp_object for hierarchical beta
library(dirichletprocess)
source("tests/testthat/helper-testing.R")

# Create test data like in the test
group1_data <- rbeta(25, 2, 5)
group2_data <- rbeta(25, 5, 2) 
group3_data <- rbeta(25, 1, 1)
hierarchical_data <- list(group1_data, group2_data, group3_data)
attr(hierarchical_data, "maxY") <- 1

cat("Testing create_dp_object with hierarchical_beta...\n")

tryCatch({
  dp <- create_dp_object("hierarchical_beta", hierarchical_data)
  cat("create_dp_object succeeded\n")
  cat("DP class:", class(dp), "\n")
  cat("Number of individual DPs:", length(dp$indDP), "\n")
  
  # Check individual DP structure
  cat("Individual DP 1 class:", class(dp$indDP[[1]]), "\n")
  cat("Individual DP 1 has aux:", "aux" %in% names(dp$indDP[[1]]), "\n")
  
  # Try one MCMC step
  cat("Attempting Fit...\n")
  dp_fitted <- Fit(dp, its = 1, updatePrior = TRUE)
  cat("Fit succeeded\n")
  
}, error = function(e) {
  cat("create_dp_object failed with error:", e$message, "\n")
  print(e)
})

cat("Debugging complete\n")