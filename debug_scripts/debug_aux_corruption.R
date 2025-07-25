# Debug aux parameter corruption step by step
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

# Test each iteration separately
for (iter in 1:10) {
  cat("=== Iteration", iter, "===\n")
  
  # Check aux structure before fitting
  for (i in 1:3) {
    if ("aux" %in% names(dp_r$indDP[[i]]) && length(dp_r$indDP[[i]]$aux) > 0) {
      cat("Before: DP", i, "aux[[1]] names:", names(dp_r$indDP[[i]]$aux[[1]]), "\n")
    } else {
      cat("Before: DP", i, "aux missing or empty\n")
    }
  }
  
  tryCatch({
    dp_r <- Fit(dp_r, its = 1, updatePrior = TRUE)
    cat("Iteration", iter, "succeeded\n")
    
    # Check aux structure after fitting
    for (i in 1:3) {
      if ("aux" %in% names(dp_r$indDP[[i]]) && length(dp_r$indDP[[i]]$aux) > 0) {
        cat("After: DP", i, "aux[[1]] names:", names(dp_r$indDP[[i]]$aux[[1]]), "\n")
      } else {
        cat("After: DP", i, "aux missing or empty\n")
      }
    }
    
  }, error = function(e) {
    cat("Iteration", iter, "failed with error:", e$message, "\n")
    
    # Debug which DP and aux parameter is problematic
    for (i in 1:3) {
      cat("DP", i, "debugging:\n")
      if ("aux" %in% names(dp_r$indDP[[i]])) {
        cat("  Has aux:", length(dp_r$indDP[[i]]$aux), "items\n")
        for (j in 1:min(3, length(dp_r$indDP[[i]]$aux))) {
          if (is.null(dp_r$indDP[[i]]$aux[[j]])) {
            cat("  aux[[", j, "]] is NULL\n")
          } else if (is.list(dp_r$indDP[[i]]$aux[[j]])) {
            cat("  aux[[", j, "]] names:", names(dp_r$indDP[[i]]$aux[[j]]), "\n")
          } else {
            cat("  aux[[", j, "]] is not a list, class:", class(dp_r$indDP[[i]]$aux[[j]]), "\n")
          }
        }
      } else {
        cat("  No aux component\n")
      }
    }
    break
  })
  
  cat("\n")
}

cat("Debugging complete\n")