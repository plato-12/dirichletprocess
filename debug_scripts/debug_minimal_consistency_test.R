# Debug minimal consistency test
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

cat("Testing create_dp_object...\n")
tryCatch({
  dp_r <- create_dp_object("hierarchical_beta", hierarchical_data)
  cat("create_dp_object succeeded\n")
  
  cat("Testing single Fit iteration...\n")
  dp_r <- Fit(dp_r, its = 1, updatePrior = TRUE)
  cat("First iteration succeeded\n")
  
  cat("Testing 5 iterations...\n")
  dp_r <- Fit(dp_r, its = 5, updatePrior = TRUE)
  cat("Multiple iterations succeeded\n")
  
}, error = function(e) {
  cat("Test failed with error:", e$message, "\n")
  print(e)
  
  # Additional debugging
  cat("\nDebugging info:\n")
  if (exists("dp_r")) {
    cat("DP created successfully\n")
    cat("Number of individual DPs:", length(dp_r$indDP), "\n")
    
    for (i in 1:length(dp_r$indDP)) {
      cat("DP", i, "class:", class(dp_r$indDP[[i]]), "\n")
      cat("DP", i, "mixing distribution class:", class(dp_r$indDP[[i]]$mixingDistribution), "\n")
      if ("aux" %in% names(dp_r$indDP[[i]])) {
        cat("DP", i, "has aux parameters\n")
        if (length(dp_r$indDP[[i]]$aux) > 0) {
          cat("DP", i, "aux[[1]] names:", names(dp_r$indDP[[i]]$aux[[1]]), "\n")
        }
      } else {
        cat("DP", i, "does NOT have aux parameters\n")
      }
    }
  } else {
    cat("DP was not created\n")
  }
})

cat("Debugging complete\n")