---
output:
  pdf_document: default
  html_document: default
---
# C++ Covariance Models Usage Guide

This document demonstrates the usage of all covariance models in the Dirichlet Process package with C++ implementations enabled.

## Table of Contents

1. [Setup and Configuration](#setup-and-configuration)
2. [Univariate Models (1D)](#univariate-models-1d)
3. [Multivariate Models (2D+)](#multivariate-models-2d)
4. [Constrained Covariance Models](#constrained-covariance-models)
5. [Model Comparison Examples](#model-comparison-examples)
6. [Performance Tips](#performance-tips)
7. [Troubleshooting](#troubleshooting)

## Setup and Configuration

### Loading the Package with C++ Support

```r
library(dirichletprocess)

# Load development version (if applicable)
devtools::load_all()

# Enable C++ implementations
set_use_cpp(TRUE)

# Verify C++ is enabled
cat("C++ status:", ifelse(using_cpp(), "ENABLED", "DISABLED"), "\n")
```

### Basic Data Generation

```r
# Set seed for reproducibility
set.seed(42)

# Generate 1D data (univariate)
n1d <- 100
data_1d <- matrix(c(
  rnorm(n1d/2, mean = -2, sd = 1),
  rnorm(n1d/2, mean = 2, sd = 1)
), ncol = 1)

# Generate 2D data (bivariate)
n2d <- 150
library(mvtnorm)
data_2d <- rbind(
  rmvnorm(n2d/3, mean = c(-2, -2), sigma = diag(2)),
  rmvnorm(n2d/3, mean = c(0, 2), sigma = diag(2)),
  rmvnorm(n2d/3, mean = c(2, -1), sigma = diag(2))
)

# Generate 5D data (multivariate)
n5d <- 200
data_5d <- rbind(
  rmvnorm(n5d/2, mean = rep(-1, 5), sigma = diag(5) * 0.5),
  rmvnorm(n5d/2, mean = rep(1, 5), sigma = diag(5) * 0.5)
)
```

## Univariate Models (1D)

### E Model (Equal Variance)

The E model assumes equal variances across all clusters.

```r
# Create E model prior parameters
prior_e <- list(
  mu0 = 0,          # Prior mean
  kappa0 = 1,       # Prior precision for mean
  nu = 3,           # Degrees of freedom
  Lambda = matrix(1), # Prior scale matrix
  covModel = "E"
)

# Create mixing distribution and Dirichlet process
md_e <- MvnormalCreate(prior_e)
dp_e <- DirichletProcessCreate(data_1d, md_e)

# Initialize with single cluster
dp_e <- Initialise(dp_e, numInitialClusters = 1)

# Run MCMC with C++ acceleration
dp_e_fitted <- Fit(dp_e, its = 100, progressBar = TRUE)

# Results
cat("E Model Results:\n")
cat("  Clusters found:", dp_e_fitted$numberClusters, "\n")
cat("  Final log likelihood:", tail(dp_e_fitted$logLikelihood, 1), "\n")
```

### V Model (Variable Variance)

The V model allows different variances across clusters.

```r
# Create V model prior parameters
prior_v <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1),
  covModel = "V"
)

# Create and fit V model
md_v <- MvnormalCreate(prior_v)
dp_v <- DirichletProcessCreate(data_1d, md_v)
dp_v <- Initialise(dp_v, numInitialClusters = 1)
dp_v_fitted <- Fit(dp_v, its = 100, progressBar = TRUE)

cat("V Model Results:\n")
cat("  Clusters found:", dp_v_fitted$numberClusters, "\n")
cat("  Final log likelihood:", tail(dp_v_fitted$logLikelihood, 1), "\n")
```

### FULL Model (Unrestricted - 1D)

For 1D data, FULL model is equivalent to V model but uses different parameterization.

```r
# Create FULL model for 1D data
prior_full_1d <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1),
  covModel = "FULL"
)

# Create and fit FULL model
md_full_1d <- MvnormalCreate(prior_full_1d)
dp_full_1d <- DirichletProcessCreate(data_1d, md_full_1d)
dp_full_1d <- Initialise(dp_full_1d, numInitialClusters = 1)
dp_full_1d_fitted <- Fit(dp_full_1d, its = 100, progressBar = TRUE)

cat("FULL Model (1D) Results:\n")
cat("  Clusters found:", dp_full_1d_fitted$numberClusters, "\n")
cat("  Final log likelihood:", tail(dp_full_1d_fitted$logLikelihood, 1), "\n")
```

## Multivariate Models (2D+)

### FULL Model (Unrestricted Covariance)

The FULL model places no restrictions on covariance matrices.

```r
# Create FULL model for multivariate data
d <- ncol(data_2d)  # Number of dimensions
prior_full <- list(
  mu0 = rep(0, d),
  kappa0 = 1,
  nu = d + 2,
  Lambda = diag(d),
  covModel = "FULL"
)

# Create and fit FULL model
md_full <- MvnormalCreate(prior_full)
dp_full <- DirichletProcessCreate(data_2d, md_full)
dp_full <- Initialise(dp_full, numInitialClusters = 1)
dp_full_fitted <- Fit(dp_full, its = 100, progressBar = TRUE)

cat("FULL Model (2D) Results:\n")
cat("  Clusters found:", dp_full_fitted$numberClusters, "\n")
cat("  Dimensions:", d, "\n")
cat("  Final log likelihood:", tail(dp_full_fitted$logLikelihood, 1), "\n")

# Examine cluster parameters
for (i in 1:dp_full_fitted$numberClusters) {
  cat("Cluster", i, "covariance matrix:\n")
  print(dp_full_fitted$clusterParameters[[i]]$sig[,,i])
}
```

## Constrained Covariance Models

### EII Model (Equal, Isotropic, Identical)

All clusters have the same spherical covariance matrix.

```r
# Create EII model
prior_eii <- list(
  mu0 = rep(0, ncol(data_2d)),
  kappa0 = 1,
  nu = ncol(data_2d) + 2,
  Lambda = diag(ncol(data_2d)),
  covModel = "EII"
)

md_eii <- MvnormalCreate(prior_eii)
dp_eii <- DirichletProcessCreate(data_2d, md_eii)
dp_eii <- Initialise(dp_eii, numInitialClusters = 1)
dp_eii_fitted <- Fit(dp_eii, its = 100, progressBar = TRUE)

cat("EII Model Results:\n")
cat("  Clusters found:", dp_eii_fitted$numberClusters, "\n")
cat("  Model type: Equal, Isotropic, Identical orientation\n")
```

### VII Model (Variable Volume, Isotropic, Identical)

Clusters have different volumes but same spherical shape.

```r
# Create VII model
prior_vii <- list(
  mu0 = rep(0, ncol(data_2d)),
  kappa0 = 1,
  nu = ncol(data_2d) + 2,
  Lambda = diag(ncol(data_2d)),
  covModel = "VII"
)

md_vii <- MvnormalCreate(prior_vii)
dp_vii <- DirichletProcessCreate(data_2d, md_vii)
dp_vii <- Initialise(dp_vii, numInitialClusters = 1)
dp_vii_fitted <- Fit(dp_vii, its = 100, progressBar = TRUE)

cat("VII Model Results:\n")
cat("  Clusters found:", dp_vii_fitted$numberClusters, "\n")
cat("  Model type: Variable volume, Isotropic, Identical orientation\n")
```

### EEI Model (Equal Volume, Equal Shape, Isotropic)

```r
# Create EEI model
prior_eei <- list(
  mu0 = rep(0, ncol(data_2d)),
  kappa0 = 1,
  nu = ncol(data_2d) + 2,
  Lambda = diag(ncol(data_2d)),
  covModel = "EEI"
)

md_eei <- MvnormalCreate(prior_eei)
dp_eei <- DirichletProcessCreate(data_2d, md_eei)
dp_eei <- Initialise(dp_eei, numInitialClusters = 1)
dp_eei_fitted <- Fit(dp_eei, its = 100, progressBar = TRUE)

cat("EEI Model Results:\n")
cat("  Clusters found:", dp_eei_fitted$numberClusters, "\n")
```

### VEI Model (Variable Volume, Equal Shape, Isotropic)

```r
# Create VEI model
prior_vei <- list(
  mu0 = rep(0, ncol(data_2d)),
  kappa0 = 1,
  nu = ncol(data_2d) + 2,
  Lambda = diag(ncol(data_2d)),
  covModel = "VEI"
)

md_vei <- MvnormalCreate(prior_vei)
dp_vei <- DirichletProcessCreate(data_2d, md_vei)
dp_vei <- Initialise(dp_vei, numInitialClusters = 1)
dp_vei_fitted <- Fit(dp_vei, its = 100, progressBar = TRUE)

cat("VEI Model Results:\n")
cat("  Clusters found:", dp_vei_fitted$numberClusters, "\n")
```

### EVI Model (Equal Volume, Variable Shape, Isotropic)

```r
# Create EVI model
prior_evi <- list(
  mu0 = rep(0, ncol(data_2d)),
  kappa0 = 1,
  nu = ncol(data_2d) + 2,
  Lambda = diag(ncol(data_2d)),
  covModel = "EVI"
)

md_evi <- MvnormalCreate(prior_evi)
dp_evi <- DirichletProcessCreate(data_2d, md_evi)
dp_evi <- Initialise(dp_evi, numInitialClusters = 1)
dp_evi_fitted <- Fit(dp_evi, its = 100, progressBar = TRUE)

cat("EVI Model Results:\n")
cat("  Clusters found:", dp_evi_fitted$numberClusters, "\n")
```

### VVI Model (Variable Volume, Variable Shape, Isotropic)

```r
# Create VVI model
prior_vvi <- list(
  mu0 = rep(0, ncol(data_2d)),
  kappa0 = 1,
  nu = ncol(data_2d) + 2,
  Lambda = diag(ncol(data_2d)),
  covModel = "VVI"
)

md_vvi <- MvnormalCreate(prior_vvi)
dp_vvi <- DirichletProcessCreate(data_2d, md_vvi)
dp_vvi <- Initialise(dp_vvi, numInitialClusters = 1)
dp_vvi_fitted <- Fit(dp_vvi, its = 100, progressBar = TRUE)

cat("VVI Model Results:\n")
cat("  Clusters found:", dp_vvi_fitted$numberClusters, "\n")
```

## Model Comparison Examples

### Comparing Models on Same Dataset

```r
# Function to fit and evaluate model
fit_and_evaluate <- function(data, prior_params, model_name) {
  md <- MvnormalCreate(prior_params)
  dp <- DirichletProcessCreate(data, md)
  dp <- Initialise(dp, numInitialClusters = 1)
  
  start_time <- Sys.time()
  dp_fitted <- Fit(dp, its = 50, progressBar = FALSE)
  end_time <- Sys.time()
  
  return(list(
    model = model_name,
    clusters = dp_fitted$numberClusters,
    log_likelihood = tail(dp_fitted$logLikelihood, 1),
    time_sec = as.numeric(end_time - start_time),
    fitted_object = dp_fitted
  ))
}

# Compare models on 2D data
models_2d <- list(
  "FULL" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "FULL"),
  "EII" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "EII"),
  "VII" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "VII"),
  "EEI" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "EEI"),
  "VEI" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "VEI"),
  "EVI" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "EVI"),
  "VVI" = list(mu0 = rep(0, 2), kappa0 = 1, nu = 4, Lambda = diag(2), covModel = "VVI")
)

# Fit all models
cat("Comparing all 2D covariance models...\n")
results_2d <- list()
for (model_name in names(models_2d)) {
  cat("Fitting", model_name, "model...")
  results_2d[[model_name]] <- fit_and_evaluate(data_2d, models_2d[[model_name]], model_name)
  cat(" ✓ Done\n")
}

# Display comparison results
cat("\n=== MODEL COMPARISON RESULTS ===\n")
comparison_df <- data.frame(
  Model = names(results_2d),
  Clusters = sapply(results_2d, function(x) x$clusters),
  LogLikelihood = sapply(results_2d, function(x) round(x$log_likelihood, 2)),
  Time_Sec = sapply(results_2d, function(x) round(x$time_sec, 2))
)

# Sort by log likelihood (higher is better)
comparison_df <- comparison_df[order(-comparison_df$LogLikelihood), ]
print(comparison_df)

cat("\nBest model by log likelihood:", comparison_df$Model[1], "\n")
cat("Fastest model:", comparison_df$Model[which.min(comparison_df$Time_Sec)], "\n")
```

### High-Dimensional Example (5D)

```r
# Test models on 5D data
d5 <- 5
prior_5d <- list(
  mu0 = rep(0, d5),
  kappa0 = 1,
  nu = d5 + 2,
  Lambda = diag(d5)
)

# Test FULL and constrained models
models_5d <- c("FULL", "EII", "VII", "VVI")

cat("Testing models on 5D data...\n")
for (model in models_5d) {
  cat("Testing", model, "model...")
  
  prior_5d$covModel <- model
  result <- fit_and_evaluate(data_5d, prior_5d, model)
  
  cat(sprintf(" ✓ %d clusters, %.2f sec\n", result$clusters, result$time_sec))
}
```

## Performance Tips

### 1. Optimal MCMC Settings

```r
# For development/testing
dp_fitted <- Fit(dp, its = 50, progressBar = TRUE)

# For research quality
dp_fitted <- Fit(dp, its = 200, progressBar = TRUE)

# For production/final results
dp_fitted <- Fit(dp, its = 1000, progressBar = TRUE)
```

### 2. Model Selection Guidelines

```r
# For 1D data
# - Use E model for equal variance assumption
# - Use V model for variable variance
# - Use FULL for maximum flexibility

# For 2D+ data
# - Start with EII for simple spherical clusters
# - Use VII if clusters have different sizes
# - Use FULL for maximum flexibility (slower)
# - Use VVI for moderate flexibility

# Rule of thumb:
# - Small datasets (n < 100): Use constrained models (EII, VII, etc.)
# - Large datasets (n > 500): FULL model becomes more reliable
# - High dimensions (d > 10): Prefer constrained models
```

### 3. Prior Parameter Guidelines

```r
# Conservative priors (less informative)
prior_conservative <- list(
  mu0 = rep(0, d),      # Center at origin
  kappa0 = 0.1,         # Low precision (more uncertainty)
  nu = d + 1,           # Minimal degrees of freedom
  Lambda = diag(d) * 0.1, # Small scale
  covModel = "FULL"
)

# Informative priors (when you have domain knowledge)
prior_informative <- list(
  mu0 = c(1, -1),       # Known approximate cluster centers
  kappa0 = 2,           # Higher precision
  nu = d + 5,           # More degrees of freedom
  Lambda = diag(d) * 2, # Larger scale
  covModel = "FULL"
)
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Model Compatibility

```r
# ✗ WRONG: E/V models with multivariate data
# prior_wrong <- list(covModel = "E", mu0 = rep(0, 3), ...)  # Will fail

# ✓ CORRECT: Check dimensions first
check_model_compatibility <- function(data, cov_model) {
  d <- ncol(data)
  
  if (cov_model %in% c("E", "V") && d > 1) {
    stop("Models E and V are only for univariate data (d=1)")
  }
  
  if (cov_model %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI") && d == 1) {
    stop(paste("Model", cov_model, "requires multivariate data (d>1)"))
  }
  
  return(TRUE)
}

# Use before creating models
check_model_compatibility(data_2d, "EII")  # ✓ OK
# check_model_compatibility(data_1d, "EII")  # ✗ Would error
```

#### 2. C++ Availability Check

```r
# Check if C++ is working properly
diagnose_cpp <- function() {
  cat("C++ Status Check:\n")
  cat("  using_cpp():", using_cpp(), "\n")
  cat("  can_use_cpp():", can_use_cpp(), "\n")
  
  # Test simple model
  tryCatch({
    test_data <- matrix(rnorm(50), ncol = 1)
    prior <- list(mu0 = 0, kappa0 = 1, nu = 3, Lambda = matrix(1), covModel = "E")
    md <- MvnormalCreate(prior)
    dp <- DirichletProcessCreate(test_data, md)
    dp <- Initialise(dp, numInitialClusters = 1)
    result <- Fit(dp, its = 5, progressBar = FALSE)
    
    cat("  C++ test: SUCCESS\n")
    cat("  Test clusters found:", result$numberClusters, "\n")
    
  }, error = function(e) {
    cat("  C++ test: FAILED -", e$message, "\n")
    cat("  Try: set_use_cpp(FALSE) to use R implementation\n")
  })
}

# Run diagnosis
diagnose_cpp()
```

#### 3. Memory and Performance

```r
# For large datasets, use batch processing
fit_large_dataset <- function(data, prior_params, batch_size = 1000) {
  if (nrow(data) <= batch_size) {
    # Small dataset - process normally
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(data, md)
    dp <- Initialise(dp, numInitialClusters = 1)
    return(Fit(dp, its = 100, progressBar = TRUE))
  } else {
    # Large dataset - could implement subsampling or other strategies
    cat("Large dataset detected. Consider subsampling or using fewer MCMC iterations.\n")
    sample_idx <- sample(nrow(data), batch_size)
    data_subset <- data[sample_idx, , drop = FALSE]
    
    md <- MvnormalCreate(prior_params)
    dp <- DirichletProcessCreate(data_subset, md)
    dp <- Initialise(dp, numInitialClusters = 1)
    return(Fit(dp, its = 100, progressBar = TRUE))
  }
}
```

## Summary

This guide demonstrates the usage of all covariance models available in the Dirichlet Process package with C++ acceleration:

- **Univariate models (1D)**: E, V, FULL
- **Multivariate models (2D+)**: FULL, EII, VII, EEI, VEI, EVI, VVI

Key points:
1. Always enable C++ with `set_use_cpp(TRUE)` for better performance
2. Choose appropriate models based on data dimensionality
3. Use constrained models for stability with smaller datasets
4. Use FULL model for maximum flexibility with larger datasets
5. Monitor convergence through log likelihood and cluster counts

For more information, see the package documentation and vignettes.
