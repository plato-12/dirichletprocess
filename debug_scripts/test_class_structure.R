# Test class structure
library(dirichletprocess)
set_use_cpp(FALSE)  # Turn off C++

cat("=== TESTING CLASS STRUCTURE ===\n")

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

cat("Classes of dp object:\n")
print(class(dp))

cat("\nMethods available for Initialise:\n")
print(methods("Initialise"))

# Check which method will be called
cat("\nWhich method will be called?\n")
tryCatch({
  method <- getS3method("Initialise", class(dp))
  cat(sprintf("Method: %s\n", deparse(method)))
}, error = function(e) {
  cat(sprintf("Method lookup failed: %s\n", e$message))
})

# Test with class structure
cat("\nTesting inherits function:\n")
cat(sprintf("inherits(dp, 'mvnormal'): %s\n", inherits(dp, 'mvnormal')))
cat(sprintf("inherits(dp, 'conjugate'): %s\n", inherits(dp, 'conjugate')))

cat("\n=== CLASS STRUCTURE TEST COMPLETE ===\n")