# Debug Predictive Probability Step by Step
library(dirichletprocess)

cat("=== Step by Step Predictive Probability Debug ===\n")

# Test setup
set.seed(123)
x <- -0.5604756  # First test data point
prior_params <- c(0, 1, 1, 1)  # mu0, kappa0, alpha0, beta0

mu0 <- prior_params[1]
kappa0 <- prior_params[2]
alpha0 <- prior_params[3]
beta0 <- prior_params[4]

cat("Input x:", x, "\n")
cat("Prior parameters:", prior_params, "\n")

# Manual calculation of posterior parameters for single observation
n.x <- 1
ybar <- x

mu.n <- (kappa0 * mu0 + n.x * ybar)/(kappa0 + n.x)
kappa.n <- kappa0 + n.x
alpha.n <- alpha0 + n.x/2
beta.n <- beta0 + 0.5 * sum((x - ybar)^2) + kappa0 * n.x * (ybar - mu0)^2/(2 * (kappa0 + n.x))

cat("Manual posterior parameters:\n")
cat("  mu.n:", mu.n, "\n")
cat("  kappa.n:", kappa.n, "\n")
cat("  alpha.n:", alpha.n, "\n")
cat("  beta.n:", beta.n, "\n")

# Manual predictive calculation using R's formula
manual_predictive <- (gamma(alpha.n)/gamma(alpha0)) *
  ((beta0^alpha0)/(beta.n^alpha.n)) *
  sqrt(kappa0/kappa.n)

cat("Manual predictive probability:", manual_predictive, "\n")

# Now check what R's PosteriorParameters function returns
dp <- DirichletProcessGaussian(c(0))  # Dummy DP with same priors
posterior_params_r <- PosteriorParameters(dp$mixingDistribution, x)
cat("R PosteriorParameters result:\n")
cat("  [1]:", posterior_params_r[1], "\n")
cat("  [2]:", posterior_params_r[2], "\n")
cat("  [3]:", posterior_params_r[3], "\n")
cat("  [4]:", posterior_params_r[4], "\n")

# R predictive calculation
r_predictive <- (gamma(posterior_params_r[3])/gamma(alpha0)) *
  ((beta0^alpha0)/(posterior_params_r[4]^posterior_params_r[3])) *
  sqrt(kappa0/posterior_params_r[2])

cat("R-based predictive probability:", r_predictive, "\n")

# Check the actual Predictive function result
predictive_result <- Predictive(dp$mixingDistribution, x)
cat("R Predictive function result:", predictive_result, "\n")

# Compare all results
cat("\n=== Comparison ===\n")
cat("Manual calculation:", manual_predictive, "\n")
cat("R-based calculation:", r_predictive, "\n")
cat("R Predictive function:", predictive_result, "\n")
cat("All match:", abs(manual_predictive - r_predictive) < 1e-10 && abs(r_predictive - predictive_result) < 1e-10, "\n")