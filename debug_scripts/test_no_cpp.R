# Test without C++
library(dirichletprocess)
set_use_cpp(FALSE)  # Turn off C++

cat("=== TESTING WITHOUT C++ ===\n")

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
dp <- DirichletProcessCreate(data_1d, md)

cat("Calling Initialise with C++ disabled...\n")
tryCatch({
  dp1 <- Initialise(dp, numInitialClusters = 2)
  cat("SUCCESS without C++\n")
}, error = function(e) {
  cat(sprintf("FAILED without C++: %s\n", e$message))
})

cat("\n=== NO C++ TEST COMPLETE ===\n")