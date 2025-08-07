
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dirichletprocess

[![R build
status](https://github.com/dm13450/dirichletprocess/workflows/R-CMD-check/badge.svg)](https://github.com/dm13450/dirichletprocess/actions)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/dm13450/dirichletprocess?branch=master&svg=true)](https://ci.appveyor.com/project/dm13450/dirichletprocess)
[![Coverage
Status](https://codecov.io/gh/dm13450/dirichletprocess/branch/master/graph/badge.svg)](https://app.codecov.io/gh/dm13450/dirichletprocess)

The dirichletprocess package provides tools for building custom
Dirichlet process mixture models for nonparametric Bayesian analysis. 
The package features high-performance C++ implementations alongside pure R 
implementations, offering significant speed improvements while maintaining 
full compatibility and automatic fallback mechanisms.

**Key Features:**
- Pre-built distributions: Normal, Beta, Exponential, Weibull, Multivariate Normal
- Hierarchical Dirichlet Process models
- High-performance C++ backend with automatic R fallback
- Comprehensive covariance model support (EII, VII, EEI, VEI, EVI, VVI, FULL)
- Advanced MCMC algorithms (Neal's Algorithm 4 & 8)
- Extensive validation and testing framework

Model your data nonparametrically in as little as four lines of code.

## Installation

You can install the stable release of dirichletprocess from CRAN:

``` r
install.packages("dirichletprocess")
```

You can also install the development build of dirichletprocess from
github with:

``` r
# install.packages("devtools")
devtools::install_github("dm13450/dirichletprocess")
```

For a full guide to the package and its capabilities please consult the
vignette:

``` r
browseVignettes(package = "dirichletprocess")
```

## Examples

### Density Estimation

Dirichlet processes can be used for nonparametric density estimation.

``` r
faithfulTransformed <- faithful$waiting - mean(faithful$waiting)
faithfulTransformed <- faithfulTransformed/sd(faithful$waiting)

# Using R implementation (default)
dp <- DirichletProcessGaussian(faithfulTransformed, cpp = FALSE)
dp <- Fit(dp, 100, progressBar = FALSE)
plot(dp)

# For better performance, use C++ implementation
# dp <- DirichletProcessGaussian(faithfulTransformed, cpp = TRUE)
```

<img src=https://github.com/dm13450/dirichletprocess/raw/master/vignettes/img/density-1.png width=50% />

### Clustering

Dirichlet processes can also be used to cluster data based on their
common distribution parameters.

``` r
faithfulTrans <- scale(faithful)

# Using R implementation (default)
dpCluster <- DirichletProcessMvnormal(faithfulTrans, cpp = FALSE)
dpCluster <- Fit(dpCluster, 2000, progressBar = FALSE)
plot(dpCluster)

# For better performance with large datasets, use C++ implementation
# dpCluster <- DirichletProcessMvnormal(faithfulTrans, cpp = TRUE)
```

<img src=https://github.com/dm13450/dirichletprocess/raw/master/vignettes/img/clustering-1.png width=50% />

For more detailed explanations and examples see the vignette.

## Performance & Implementation Control

The package provides both R and high-performance C++ implementations for all distributions. You can now explicitly control which implementation to use with the `cpp` parameter in all constructor functions.

### Using the cpp Parameter

**New Feature (v0.5.0+):** All distribution constructors now accept a `cpp` parameter to explicitly choose the implementation:

```r
library(dirichletprocess)
y <- rt(200, 3) + 2

# Use R implementation (default)
dp_r <- DirichletProcessGaussian(y, cpp = FALSE)
dp_r <- Fit(dp_r, 1000)

# Use C++ implementation for better performance
dp_cpp <- DirichletProcessGaussian(y, cpp = TRUE)
dp_cpp <- Fit(dp_cpp, 1000)
```

### Available for All Distributions

The `cpp` parameter works with **all** distribution types:

```r
# Normal distributions
dp <- DirichletProcessGaussian(data, cpp = TRUE)
dp <- DirichletProcessGaussianFixedVariance(data, sigma = 1, cpp = TRUE)

# Other distributions
dp <- DirichletProcessBeta(data, cpp = TRUE)
dp <- DirichletProcessExponential(data, cpp = TRUE)
dp <- DirichletProcessWeibull(data, g0Priors = c(1, 1, 1, 1), cpp = TRUE)

# Multivariate distributions
dp <- DirichletProcessMvnormal(mvdata, cpp = TRUE)
dp <- DirichletProcessMvnormal2(mvdata, cpp = TRUE)

# Hierarchical models
dp <- DirichletProcessHierarchicalBeta(dataList, maxY = 1, cpp = TRUE)
dp <- DirichletProcessHierarchicalMvnormal2(dataList, cpp = TRUE)

# Markov models
dp <- DirichletHMMCreate(data, mdobj, alpha = 1, beta = 1, cpp = TRUE)
```

### Global Control (Legacy Method)

You can still control the implementation globally using the legacy functions:

```r
# Check current implementation preference
using_cpp()  # Returns TRUE if C++ backend is preferred

# Set global preference (affects all new objects)
set_use_cpp(TRUE)   # Prefer C++ implementations
set_use_cpp(FALSE)  # Prefer R implementations
```

### Performance Benefits

C++ implementations provide substantial speedups for large datasets while maintaining **identical results** to R implementations:

- **~2-10x faster** for most distributions
- **Automatic fallback** to R if C++ unavailable  
- **Identical statistical results** guaranteed
- **Memory efficient** for large datasets

### Default Behavior

- **Default**: `cpp = FALSE` (R implementation) for predictable, cross-platform behavior
- **Recommendation**: Use `cpp = TRUE` for large datasets or production workflows
- **Compatibility**: Both implementations produce identical statistical results

## Supported Distributions

**Conjugate Models:**
- Normal (Gaussian) with Inverse-Gamma prior
- Exponential with Gamma prior  
- Multivariate Normal with Normal-Wishart prior (all covariance models)

**Non-Conjugate Models:**
- Beta with Uniform priors
- Weibull with Uniform priors
- Multivariate Normal with semi-conjugate priors

**Hierarchical Models:**
- Hierarchical Beta
- Hierarchical Multivariate Normal (two variants)

## Covariance Models

For multivariate normal distributions, the package supports:
- **FULL**: Unrestricted covariance matrices
- **EII, VII, EEI, VEI, EVI, VVI**: Constrained covariance models

### Tutorials

I've written a number of tutorials:

-   [Non parametric
    priors](https://dm13450.github.io/2019/02/22/Nonparametric-Prior.html)
-   [Calculating cluster
    probabilities](https://dm13450.github.io/2018/11/21/Cluster-Probabilities.html)
-   [Clustering](https://dm13450.github.io/2018/05/30/Clustering.html)
-   [Custom
    mixtures](https://dm13450.github.io/2018/02/21/Custom-Distributions-Conjugate.html)
-   [Density
    estimation](https://dm13450.github.io/2018/02/01/Dirichlet-Density.html)
-   [Checking
    convergence](https://dm13450.github.io/2020/01/11/Dirichlet-Convergence.html)

and some case studies:

-   [State of the Market - Infinite State Hidden Markov
    Models](https://dm13450.github.io/2020/06/03/State-of-the-Market.html)
-   [Palmer Penguins and an Introduction to Dirichlet
    Processes](https://dm13450.github.io/2020/09/28/PriorToPosterior.html)
