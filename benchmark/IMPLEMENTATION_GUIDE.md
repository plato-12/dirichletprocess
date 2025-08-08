# Implementation Guide: R vs C++ in dirichletprocess

This guide explains how to control and navigate between R and C++ implementations in the dirichletprocess package.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Understanding Implementations](#understanding-implementations)
3. [Per-Object Control (Recommended)](#per-object-control-recommended)
4. [Global Control (Legacy)](#global-control-legacy)
5. [Checking Current Implementation](#checking-current-implementation)
6. [All Supported Distributions](#all-supported-distributions)
7. [Performance Considerations](#performance-considerations)
8. [Troubleshooting](#troubleshooting)
9. [Migration from Older Versions](#migration-from-older-versions)

## Quick Start

### Basic Usage with Explicit Control (New Method)

```r
library(dirichletprocess)
y <- rt(200, 3) + 2

# Use R implementation (default, stable)
dp_r <- DirichletProcessGaussian(y, cpp = FALSE)
dp_r <- Fit(dp_r, 1000)

# Use C++ implementation (faster)
dp_cpp <- DirichletProcessGaussian(y, cpp = TRUE)
dp_cpp <- Fit(dp_cpp, 1000)
```

### Check What Implementation You're Using

```r
# Check global preference
using_cpp()

# This tells you if the package will prefer C++ or R for new objects
# TRUE = C++ preferred, FALSE = R preferred
```

## Understanding Implementations

### R Implementation
- **Pros**: 
  - More stable across platforms
  - Easier to debug
  - Works on all systems
  - Default behavior (predictable)
- **Cons**: 
  - Slower for large datasets
  - Higher memory usage

### C++ Implementation  
- **Pros**:
  - 2-10x faster performance
  - Memory efficient
  - Identical statistical results
- **Cons**:
  - Requires successful C++ compilation
  - Platform dependent
  - Harder to debug if issues arise

## Per-Object Control (Recommended)

### The `cpp` Parameter

All distribution constructors now accept a `cpp` parameter:

```r
# Syntax: DistributionFunction(data, ..., cpp = TRUE/FALSE)
dp <- DirichletProcessGaussian(data, cpp = FALSE)  # Use R
dp <- DirichletProcessGaussian(data, cpp = TRUE)   # Use C++
```

### Examples for All Distribution Types

#### Basic Distributions

```r
# Normal distributions
dp1 <- DirichletProcessGaussian(data, cpp = TRUE)
dp2 <- DirichletProcessGaussianFixedVariance(data, sigma = 1, cpp = TRUE)

# Other univariate distributions
dp3 <- DirichletProcessBeta(data, cpp = TRUE)
dp4 <- DirichletProcessBeta2(data, maxY = 1, cpp = TRUE)
dp5 <- DirichletProcessExponential(data, cpp = TRUE)
dp6 <- DirichletProcessWeibull(data, g0Priors = c(1,1,1,1), cpp = TRUE)
```

#### Multivariate Distributions

```r
# Multivariate normal data
mvdata <- matrix(rnorm(200), ncol = 2)

dp7 <- DirichletProcessMvnormal(mvdata, cpp = TRUE)
dp8 <- DirichletProcessMvnormal2(mvdata, cpp = TRUE)
```

#### Hierarchical Models

```r
# Hierarchical models with multiple datasets
y1 <- runif(30, 0, 0.5)
y2 <- runif(30, 0.5, 1)
dataList <- list(y1, y2)

dp9 <- DirichletProcessHierarchicalBeta(dataList, maxY = 1, cpp = TRUE)
dp10 <- DirichletProcessHierarchicalMvnormal2(dataList, cpp = TRUE)
```

#### Markov Models

```r
# Hidden Markov Models
mdobj <- GaussianMixtureCreate(c(0, 1, 1, 1))
dp11 <- DirichletHMMCreate(data, mdobj, alpha = 1, beta = 1, cpp = TRUE)
```

### Why Use Per-Object Control?

```r
# Different strategies for different analyses
dp_exploratory <- DirichletProcessGaussian(small_data, cpp = FALSE)  # R for debugging
dp_production <- DirichletProcessGaussian(large_data, cpp = TRUE)    # C++ for speed

# Mixed workflow
dp_prototype <- DirichletProcessBeta(subset_data, cpp = FALSE)       # Test with R
dp_final <- DirichletProcessBeta(full_dataset, cpp = TRUE)           # Deploy with C++
```

## Global Control (Legacy)

### Setting Global Preferences

```r
# Set global preference for C++
set_use_cpp(TRUE)
dp1 <- DirichletProcessGaussian(data)  # Will use C++ if available

# Set global preference for R
set_use_cpp(FALSE)  
dp2 <- DirichletProcessGaussian(data)  # Will use R implementation

# Check current global setting
current_setting <- using_cpp()
print(paste("Currently using C++:", current_setting))
```

### Package Default Behavior

```r
# When you load the package, the default is:
library(dirichletprocess)
using_cpp()  # Returns FALSE (R implementation by default)
```

## Checking Current Implementation

### Available Check Functions

```r
# Check global preference
using_cpp()  # TRUE = C++ preferred, FALSE = R preferred

# Check if C++ implementations are available on your system
library(dirichletprocess)
y <- rnorm(100)
dp <- DirichletProcessGaussian(y, cpp = TRUE)

# The package will automatically fall back to R if C++ isn't available
# No error will occur - it's seamless
```

### Diagnostic Commands

```r
# Check if your system supports C++ implementations
tryCatch({
  y <- rnorm(50)
  dp <- DirichletProcessGaussian(y, cpp = TRUE)
  dp <- Fit(dp, 10, progressBar = FALSE)
  cat("C++ implementation is working!\n")
}, error = function(e) {
  cat("C++ implementation not available, using R fallback\n")
  cat("Error:", e$message, "\n")
})
```

### During Analysis

```r
# Before fitting - check your preference
cat("Global preference for C++:", using_cpp(), "\n")

# Fit your model
dp <- DirichletProcessGaussian(data, cpp = TRUE)
dp <- Fit(dp, 1000)

# After fitting - you can verify by testing performance
system.time({
  dp_r <- DirichletProcessGaussian(data, cpp = FALSE)
  dp_r <- Fit(dp_r, 100, progressBar = FALSE)
})

system.time({
  dp_cpp <- DirichletProcessGaussian(data, cpp = TRUE)  
  dp_cpp <- Fit(dp_cpp, 100, progressBar = FALSE)
})
```

## All Supported Distributions

### Complete List with cpp Parameter

All these functions support the `cpp = TRUE/FALSE` parameter:

| Distribution Type | Function | Usage |
|------------------|----------|--------|
| Normal | `DirichletProcessGaussian()` | `DirichletProcessGaussian(y, cpp = TRUE)` |
| Normal Fixed Var | `DirichletProcessGaussianFixedVariance()` | `DirichletProcessGaussianFixedVariance(y, sigma, cpp = TRUE)` |
| Beta | `DirichletProcessBeta()` | `DirichletProcessBeta(y, cpp = TRUE)` |
| Beta2 | `DirichletProcessBeta2()` | `DirichletProcessBeta2(y, maxY, cpp = TRUE)` |
| Exponential | `DirichletProcessExponential()` | `DirichletProcessExponential(y, cpp = TRUE)` |
| Weibull | `DirichletProcessWeibull()` | `DirichletProcessWeibull(y, g0Priors, cpp = TRUE)` |
| Multivariate Normal | `DirichletProcessMvnormal()` | `DirichletProcessMvnormal(y, cpp = TRUE)` |
| Multivariate Normal2 | `DirichletProcessMvnormal2()` | `DirichletProcessMvnormal2(y, cpp = TRUE)` |
| Hierarchical Beta | `DirichletProcessHierarchicalBeta()` | `DirichletProcessHierarchicalBeta(dataList, maxY, cpp = TRUE)` |
| Hierarchical MVN2 | `DirichletProcessHierarchicalMvnormal2()` | `DirichletProcessHierarchicalMvnormal2(dataList, cpp = TRUE)` |
| Markov HMM | `DirichletHMMCreate()` | `DirichletHMMCreate(x, mdobj, alpha, beta, cpp = TRUE)` |

**Total: 11 distributions with full cpp parameter support**

## Performance Considerations

### When to Use C++ (`cpp = TRUE`)

- **Large datasets** (n > 1000)
- **Production environments**
- **Repeated analyses**
- **Time-critical applications**
- **Memory-constrained systems**

```r
# Example: Large dataset analysis
large_data <- rt(10000, 3) + 2
dp_fast <- DirichletProcessGaussian(large_data, cpp = TRUE)
system.time(dp_fast <- Fit(dp_fast, 1000))  # Much faster
```

### When to Use R (`cpp = FALSE`)

- **Debugging and development**
- **Small datasets** (n < 100)  
- **Cross-platform compatibility concerns**
- **Learning and exploration**
- **Systems without C++ compilation**

```r
# Example: Exploratory analysis
small_sample <- rt(50, 3) + 2
dp_explore <- DirichletProcessGaussian(small_sample, cpp = FALSE)
dp_explore <- Fit(dp_explore, 100)
plot(dp_explore)  # Quick exploration with R
```

### Performance Comparison

```r
# Benchmark both implementations
library(microbenchmark)

data <- rt(1000, 3) + 2

benchmark <- microbenchmark(
  R_impl = {
    dp <- DirichletProcessGaussian(data, cpp = FALSE)
    Fit(dp, 100, progressBar = FALSE)
  },
  CPP_impl = {
    dp <- DirichletProcessGaussian(data, cpp = TRUE)  
    Fit(dp, 100, progressBar = FALSE)
  },
  times = 5
)

print(benchmark)
```

## Troubleshooting

### Common Issues and Solutions

#### C++ Implementation Not Available

```r
# If you get errors with cpp = TRUE, fall back to R
tryCatch({
  dp <- DirichletProcessGaussian(data, cpp = TRUE)
}, error = function(e) {
  cat("C++ not available, using R implementation\n")
  dp <- DirichletProcessGaussian(data, cpp = FALSE)
})
```

#### Checking C++ Compilation Status

```r
# Test if C++ components are properly compiled
library(dirichletprocess)

# This should work without errors if C++ is available
test_data <- rnorm(20)
test_dp <- DirichletProcessGaussian(test_data, cpp = TRUE)
test_result <- Fit(test_dp, 5, progressBar = FALSE)

if (inherits(test_result, "dirichletprocess")) {
  cat("C++ implementation is working correctly!\n")
} else {
  cat("Issue with C++ implementation\n")
}
```

#### Mixed Behavior Debugging

```r
# If you're getting unexpected behavior, check:
cat("Global setting:", using_cpp(), "\n")

# Create objects with explicit settings
dp1 <- DirichletProcessGaussian(data, cpp = TRUE)
dp2 <- DirichletProcessGaussian(data, cpp = FALSE)

# Both should work (though dp1 might fall back to R if C++ unavailable)
```

### Performance Issues

```r
# If C++ seems slower than expected:
# 1. Check data size (C++ overhead for small datasets)
# 2. Verify C++ compilation was successful
# 3. Try with larger datasets

small_data <- rnorm(50)
large_data <- rnorm(5000)

# C++ might be slower for small data due to overhead
system.time(Fit(DirichletProcessGaussian(small_data, cpp = TRUE), 100))
system.time(Fit(DirichletProcessGaussian(small_data, cpp = FALSE), 100))

# C++ should be faster for large data  
system.time(Fit(DirichletProcessGaussian(large_data, cpp = TRUE), 100))
system.time(Fit(DirichletProcessGaussian(large_data, cpp = FALSE), 100))
```

## Migration from Older Versions

### If You Previously Used Global Control

**Old approach (still works):**
```r
# Legacy method - still functional
set_use_cpp(TRUE)  
dp <- DirichletProcessGaussian(data)  # Uses global setting
```

**New recommended approach:**
```r
# Explicit per-object control (recommended)
dp <- DirichletProcessGaussian(data, cpp = TRUE)  # Clear and explicit
```

### Default Behavior Change

**Before version 0.5.0:**
- Package default was `cpp = TRUE` (C++ preferred)
- No `cpp` parameter in constructors

**Version 0.5.0+:**
- Package default is `cpp = FALSE` (R implementation default)  
- All constructors accept `cpp` parameter
- More predictable cross-platform behavior

### Updating Your Code

```r
# If your old code relied on automatic C++ usage:
# OLD (pre-0.5.0)
dp <- DirichletProcessGaussian(data)  # Used C++ by default

# NEW (0.5.0+) - specify explicitly
dp <- DirichletProcessGaussian(data, cpp = TRUE)  # Explicit C++ request
dp <- DirichletProcessGaussian(data, cpp = FALSE) # Explicit R request (default)
```

---

## Summary

- **Use `cpp = TRUE`** for performance with large datasets
- **Use `cpp = FALSE`** (default) for stability and debugging  
- **Use `using_cpp()`** to check global preferences
- **All 11 distributions** support the `cpp` parameter
- **Automatic fallback** ensures your code always works
- **Identical results** guaranteed between R and C++ implementations

This implementation control gives you the flexibility to choose the right tool for each analysis while maintaining the reliability and compatibility of the dirichletprocess package.