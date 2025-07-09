# Build and Test Instructions for C++ Hierarchical Beta DP Implementation

## Prerequisites

1. **R Development Tools**:
   - Windows: Install Rtools
   - macOS: Install Xcode Command Line Tools
   - Linux: Install r-base-dev

2. **Required R Packages**:
   ```r
   install.packages(c("Rcpp", "RcppArmadillo", "testthat", "devtools", "gtools"))
   ```

## Building the Package

### 1. Clean Previous Builds
```bash
# In the package directory
cd path/to/dirichletprocess

# Clean compiled objects
rm -f src/*.o src/*.so src/*.dll
rm -rf src-i386 src-x64
```

### 2. Generate Rcpp Exports
```r
# In R, from the package directory
Rcpp::compileAttributes()
```

### 3. Build and Install

#### Option A: Using devtools (Recommended)
```r
# In R
devtools::load_all()  # Load for development
devtools::install()   # Install the package
```3

#### Option B: Using R CMD
```bash
# From the parent directory
R CMD build dirichletprocess
R CMD INSTALL dirichletprocess_*.tar.gz
```

#### Option C: Using RStudio
- Open the project in RStudio
- Use Build > Install and Restart

## Testing the Implementation

### 1. Run Unit Tests
```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-hierarchical-beta-cpp.R")

# Run with more verbose output
testthat::test_local(reporter = "progress")
```

### 2. Run the Demonstration Script
```r
# Source the demonstration script
source("tests/test_hierarchical_beta_cpp.R")
```

### 3. Interactive Testing
```r
library(dirichletprocess)

# Check if C++ functions are available
get_cpp_status()

# Generate test data
set.seed(123)
dataList <- list(
  rbeta(20, 2, 5),
  rbeta(20, 5, 2),
  rbeta(20, 3, 3)
)

# Test with R implementation
dp_r <- DirichletProcessHierarchicalBeta(
  dataList = dataList,
  maxY = 1,
  priorParameters = c(2, 8),
  hyperPriorParameters = c(1, 0.125),
  gammaPriors = c(2, 4),
  alphaPriors = c(2, 4),
  mhStepSize = c(0.1, 0.1),
  numSticks = 20,
  mhDraws = 10
)

system.time(dp_r_fit <- Fit(dp_r, its = 100, progressBar = FALSE))

# Test with C++ implementation
enable_cpp_hierarchical_samplers(TRUE)

dp_cpp <- DirichletProcessHierarchicalBeta(
  dataList = dataList,
  maxY = 1,
  priorParameters = c(2, 8),
  hyperPriorParameters = c(1, 0.125),
  gammaPriors = c(2, 4),
  alphaPriors = c(2, 4),
  mhStepSize = c(0.1, 0.1),
  numSticks = 20,
  mhDraws = 10
)

system.time(dp_cpp_fit <- Fit(dp_cpp, its = 100, progressBar = FALSE))

# Compare performance
enable_cpp_hierarchical_samplers(FALSE)
```

## Debugging Common Issues

### 1. Compilation Errors

**Issue**: "RcppArmadillo.h not found"
```r
# Reinstall RcppArmadillo
install.packages("RcppArmadillo")
```

**Issue**: Undefined symbols
```bash
# Check the Makevars file - ensure all .cpp files are listed
cat src/Makevars
```

### 2. Runtime Errors

**Issue**: Segmentation fault
```r
# Run R with debugging
R -d gdb
# Then in gdb:
run
# After loading the package and reproducing the error:
bt  # backtrace
```

**Issue**: Index out of bounds
- Check label conversions (R uses 1-indexing, C++ uses 0-indexing)
- Verify all cluster label conversions in wrapper functions

### 3. Test Failures

**Issue**: Results don't match between R and C++
- Due to different random number generation, exact matches aren't expected
- Focus on structural similarity (same number of clusters, similar parameter ranges)

### 4. Performance Testing

```r
# Benchmark R vs C++ implementations
library(microbenchmark)

# Small dataset
small_data <- list(rbeta(10, 2, 5), rbeta(10, 5, 2))

# Create DP objects
dp_small <- DirichletProcessHierarchicalBeta(
  dataList = small_data,
  maxY = 1,
  priorParameters = c(2, 8),
  hyperPriorParameters = c(1, 0.125),
  gammaPriors = c(2, 4),
  alphaPriors = c(2, 4),
  mhStepSize = c(0.1, 0.1),
  numSticks = 10,
  mhDraws = 10
)

# Benchmark
microbenchmark(
  R = {
    enable_cpp_hierarchical_samplers(FALSE)
    Fit(dp_small, its = 10, progressBar = FALSE)
  },
  Cpp = {
    enable_cpp_hierarchical_samplers(TRUE)
    Fit(dp_small, its = 10, progressBar = FALSE)
  },
  times = 5
)
```

## Package Check

Before submitting or distributing:

```bash
# Full package check
R CMD check dirichletprocess_*.tar.gz

# Check with valgrind (Linux/macOS)
R CMD check --use-valgrind dirichletprocess_*.tar.gz

# Check for memory leaks
R -d "valgrind --leak-check=full" -e "library(dirichletprocess); source('tests/test_hierarchical_beta_cpp.R')"
```

## Continuous Integration

Add to `.travis.yml` or GitHub Actions to ensure C++ code compiles across platforms:

```yaml
# Example for GitHub Actions
- name: Check C++ compilation
  run: |
    Rscript -e "Rcpp::compileAttributes()"
    R CMD INSTALL .
    Rscript -e "devtools::test()"
```

## Notes

1. **Platform Differences**: C++ code may behave differently across platforms due to:
   - Compiler differences (gcc vs clang vs MSVC)
   - Floating-point precision
   - Random number generation

2. **Memory Management**: The C++ implementation uses RAII and smart pointers, but always verify with valgrind.

3. **Thread Safety**: Current implementation is single-threaded. For parallel MCMC, additional work would be needed.

4. **Performance**: Expect 2-10x speedup for large datasets, less improvement for small datasets due to overhead.
