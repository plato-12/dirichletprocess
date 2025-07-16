# Debug C++ function availability
library(dirichletprocess)

cat("=== C++ FUNCTION AVAILABILITY ===\n")
cat("using_cpp():", using_cpp(), "\n")
cat("using_cpp_samplers():", using_cpp_samplers(), "\n")

cat("Enabling C++ samplers...\n")
set_use_cpp(TRUE)
enable_cpp_samplers()

cat("After enabling:\n")
cat("using_cpp():", using_cpp(), "\n")
cat("using_cpp_samplers():", using_cpp_samplers(), "\n")

cat("\n=== CHECKING SPECIFIC FUNCTIONS ===\n")
functions_to_check <- c("mvnormal_likelihood_cpp", "mvnormal_posterior_draw_cpp", 
                        "mvnormal_posterior_parameters_cpp", "mvnormal_prior_draw_cpp")

for (func in functions_to_check) {
  exists_global <- exists(func)
  exists_ns <- exists(func, where = getNamespace("dirichletprocess"))
  cat(sprintf("%s: global=%s, namespace=%s\n", func, exists_global, exists_ns))
}

cat("\n=== TESTING FUNCTION CALLS ===\n")
set.seed(42)
x <- matrix(rnorm(5), nrow = 1)
mu <- rnorm(5)
sig <- diag(5)

cat("Testing mvnormal_likelihood_cpp...\n")
tryCatch({
  pkg_ns <- getNamespace("dirichletprocess")
  mvnormal_likelihood_cpp <- get("mvnormal_likelihood_cpp", pkg_ns)
  result <- mvnormal_likelihood_cpp(x, mu, sig)
  cat("Success! Result:", result, "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

cat("\n=== TESTING WRAPPER FUNCTION ===\n")
theta <- list(mu = mu, sig = sig)
priorParams <- list(covModel = "FULL")

cat("Testing mvnormal_likelihood_wrapper_cpp...\n")
tryCatch({
  pkg_ns <- getNamespace("dirichletprocess")
  mvnormal_likelihood_wrapper_cpp <- get("mvnormal_likelihood_wrapper_cpp", pkg_ns)
  result <- mvnormal_likelihood_wrapper_cpp(as.vector(x), theta, priorParams)
  cat("Success! Result:", result, "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
})