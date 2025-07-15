library(dirichletprocess)

# Get the function
mvnormal2_likelihood_cpp <- get('mvnormal2_likelihood_cpp', envir = asNamespace('dirichletprocess'))
mvnormal2_prior_draw_cpp <- get('mvnormal2_prior_draw_cpp', envir = asNamespace('dirichletprocess'))

# Create test data
priorParams <- list(
  mu0 = matrix(c(0, 0), nrow = 1),
  sigma0 = diag(2),
  phi0 = diag(2) * 2,
  nu0 = 4
)

prior_result <- mvnormal2_prior_draw_cpp(priorParams, 1)

# The C++ function expects theta as list(mu_array, sig_array) not list(mu=, sig=)
# And it expects the arrays to have the full dimensions
theta_correct <- list(prior_result$mu, prior_result$sig)

cat("Testing MVNormal2 likelihood with different input types:\n")

# Test 1: Matrix (what we were doing - should fail)
cat("1. Matrix input (should fail):\n")
test_data_matrix <- matrix(c(0.5, -0.3), nrow = 1)
tryCatch({
  result <- mvnormal2_likelihood_cpp(test_data_matrix, theta_correct)
  cat("   Success! Result:", result, "\n")
}, error = function(e) {
  cat("   Failed:", e$message, "\n")
})

# Test 2: Vector (should work)
cat("2. Vector input (should work):\n")
test_data_vector <- c(0.5, -0.3)
tryCatch({
  result <- mvnormal2_likelihood_cpp(test_data_vector, theta_correct)
  cat("   Success! Result:", result, "\n")
}, error = function(e) {
  cat("   Failed:", e$message, "\n")
})

# Test 3: As vector from matrix (should work)
cat("3. as.vector from matrix (should work):\n")
test_data_as_vector <- as.vector(test_data_matrix)
tryCatch({
  result <- mvnormal2_likelihood_cpp(test_data_as_vector, theta_correct)
  cat("   Success! Result:", result, "\n")
}, error = function(e) {
  cat("   Failed:", e$message, "\n")
})

cat("Test completed!\n")