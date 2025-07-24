# Debug Predictive Probability Calculation
library(dirichletprocess)

cat("=== Debugging Predictive Probability ===\n")

# Set up test case
set.seed(123)
test_data <- rnorm(10, mean = 0, sd = 1)
cat("Test data:", test_data, "\n")

# Create DP objects
dp_r <- DirichletProcessGaussian(test_data)
dp_cpp <- DirichletProcessGaussian(test_data)

# Check R predictive probability calculation
r_predictive <- dp_r$predictiveArray
cat("R predictive probabilities:", r_predictive, "\n")

# Get prior parameters for manual calculation
prior_params <- dp_r$mixingDistribution$priorParameters
cat("Prior parameters (mu0, kappa0, alpha0, beta0):", prior_params, "\n")

# Manual calculation of predictive probability for first data point
x1 <- test_data[1]
mu0 <- prior_params[1]
kappa0 <- prior_params[2] 
alpha0 <- prior_params[3]
beta0 <- prior_params[4]

# Predictive distribution parameters
nu <- 2 * alpha0
sigma2 <- beta0 * (kappa0 + 1) / (alpha0 * kappa0)
mu_pred <- mu0

# Student's t probability
manual_pred <- exp(lgamma((nu + 1)/2) - lgamma(nu/2) - 0.5*log(nu*pi*sigma2) - 
                  ((nu + 1)/2)*log(1 + (x1 - mu_pred)^2/(nu*sigma2)))

cat("Manual predictive calculation for x =", x1, ":", manual_pred, "\n")
cat("R predictive for same point:", r_predictive[1], "\n")
cat("Match:", abs(manual_pred - r_predictive[1]) < 1e-10, "\n")

# Check if C++ would use the same parameters
set_use_cpp(TRUE)
dp_cpp_fitted <- Fit(dp_cpp, its = 1, progressBar = FALSE)
cat("C++ fitted successfully with conjugate algorithm\n")