# Simple test for benchmark issue
library(dirichletprocess)
set_use_cpp(TRUE)

cat("=== SIMPLE BENCHMARK TEST ===\n")

# Define functions inline
generate_benchmark_data <- function(n, d, seed = 42) {
  set.seed(seed)
  if (d == 1) {
    data <- matrix(rnorm(n, mean = 0, sd = 1), ncol = 1)
  } else {
    data <- matrix(rnorm(n * d), ncol = d)
  }
  return(data)
}

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

# Test 1: Data generation
cat("Testing data generation...\n")
data_1d <- generate_benchmark_data(10, 1)
cat(sprintf("1D data: %s\n", paste(dim(data_1d), collapse = "x")))

data_2d <- generate_benchmark_data(10, 2)
cat(sprintf("2D data: %s\n", paste(dim(data_2d), collapse = "x")))

# Test 2: Parameter creation
cat("\nTesting parameter creation...\n")
params_1d <- create_prior_parameters(1, "E")
cat("1D parameters created\n")

params_2d <- create_prior_parameters(2, "FULL")
cat("2D parameters created\n")

# Test 3: Simple collection
cat("\nTesting simple collection...\n")
simple_collection <- function(model_name, data_matrix, prior_params) {
  cat(sprintf("Processing %s with %dx%d data\n", model_name, nrow(data_matrix), ncol(data_matrix)))
  
  tryCatch({
    # Create model
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(data_matrix, md)
    dp <- Initialise(dp, numInitialClusters = 2)
    
    # Quick fit
    dp <- Fit(dp, 5, progressBar = FALSE)
    
    cat(sprintf("  Success: %d clusters\n", dp$numberClusters))
    
    # Return simple metrics
    return(list(
      success = TRUE,
      n_clusters = dp$numberClusters,
      n_samples = nrow(data_matrix),
      n_features = ncol(data_matrix)
    ))
    
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(list(success = FALSE, error = e$message))
  })
}

# Test with different models
result_1d <- simple_collection("E", data_1d, params_1d)
result_2d <- simple_collection("FULL", data_2d, params_2d)

cat("\n=== RESULTS ===\n")
cat(sprintf("1D result: %s\n", result_1d$success))
cat(sprintf("2D result: %s\n", result_2d$success))

cat("\n=== SIMPLE TEST COMPLETE ===\n")