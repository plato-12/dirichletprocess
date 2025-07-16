# Test predictive issue
library(dirichletprocess)
set_use_cpp(FALSE)  # Turn off C++

cat("=== TESTING PREDICTIVE ISSUE ===\n")

# Create E model
params_1d <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1, 1, 1),
  covModel = "E"
)

md <- MvnormalCreate(params_1d)
data_1d <- matrix(rnorm(10), ncol = 1)

cat("Testing Predictive directly...\n")
tryCatch({
  pred_result <- Predictive(md, data_1d)
  cat("Predictive SUCCESS\n")
  print(pred_result)
}, error = function(e) {
  cat(sprintf("Predictive FAILED: %s\n", e$message))
})

# Test the problematic line
cat("\nTesting matrix conversion...\n")
tryCatch({
  for (i in 1:nrow(data_1d)) {
    x_i <- matrix(data_1d[i, ], nrow = 1)
    cat(sprintf("Row %d: %s\n", i, paste(x_i, collapse = ", ")))
  }
}, error = function(e) {
  cat(sprintf("Matrix conversion FAILED: %s\n", e$message))
})

cat("\n=== PREDICTIVE ISSUE TEST COMPLETE ===\n")