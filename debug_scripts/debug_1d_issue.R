# Debug 1D issue
library(dirichletprocess)
set_use_cpp(TRUE)

cat("=== DEBUGGING 1D ISSUE ===\n")

# Test parameter creation for 1D
create_prior_parameters <- function(dimensions, model_name) {
  if (model_name %in% c("E", "V")) {
    mu0 <- 0
    kappa0 <- 1
    nu <- 3
    Lambda <- matrix(1, 1, 1)
  } else {
    mu0 <- rep(0, dimensions)
    kappa0 <- 1
    nu <- dimensions + 2
    Lambda <- diag(dimensions)
  }
  
  return(list(
    mu0 = mu0,
    kappa0 = kappa0,
    nu = nu,
    Lambda = Lambda,
    covModel = model_name
  ))
}

params_1d <- create_prior_parameters(1, "E")
cat("1D parameters:\n")
print(params_1d)

# Test with manual data
data_1d <- matrix(rnorm(10), ncol = 1)
cat("\n1D data:\n")
print(str(data_1d))

# Test MvnormalCreate
cat("\nTesting MvnormalCreate...\n")
tryCatch({
  md <- MvnormalCreate(params_1d)
  cat("MvnormalCreate successful\n")
  print(md)
}, error = function(e) {
  cat(sprintf("MvnormalCreate failed: %s\n", e$message))
})

# Test DirichletProcessCreate
cat("\nTesting DirichletProcessCreate...\n")
tryCatch({
  md <- MvnormalCreate(params_1d)
  dp <- DirichletProcessCreate(data_1d, md)
  cat("DirichletProcessCreate successful\n")
  print(str(dp))
}, error = function(e) {
  cat(sprintf("DirichletProcessCreate failed: %s\n", e$message))
})

# Test Initialise
cat("\nTesting Initialise...\n")
tryCatch({
  md <- MvnormalCreate(params_1d)
  dp <- DirichletProcessCreate(data_1d, md)
  dp <- Initialise(dp, numInitialClusters = 2)
  cat("Initialise successful\n")
  print(str(dp))
}, error = function(e) {
  cat(sprintf("Initialise failed: %s\n", e$message))
})

cat("\n=== 1D DEBUG COMPLETE ===\n")