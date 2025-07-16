# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## **PROJECT MISSION: C++ IMPLEMENTATION PRIORITY**

**CRITICAL DIRECTIVE**: This project is focused on implementing and maintaining high-performance C++ implementations for the Dirichlet Process algorithms. The primary goal is to provide C++ backends that significantly improve performance over pure R implementations.

### **C++ Implementation Philosophy**

1. **Never Fall Back to R as a Solution**: When C++ implementations fail or have issues, the correct approach is to **fix the C++ code**, not to disable C++ and fall back to R
2. **Performance is Key**: C++ implementations should provide substantial performance improvements over R
3. **Correctness**: C++ implementations must produce identical results to R implementations
4. **Completeness**: All major functionality should have C++ implementations where performance matters

### **Development Guidelines for C++ Issues**

When encountering C++ implementation problems:

1. **✅ CORRECT Approach**: 
   - Analyze the root cause of the C++ issue
   - Fix the C++ code, function signatures, parameter handling, or integration
   - Test thoroughly to ensure C++ and R produce identical results
   - Document the fix and add regression tests

2. **❌ INCORRECT Approach**:
   - Disabling C++ with `set_use_cpp(FALSE)` as a solution
   - Falling back to R implementation when C++ fails
   - Using `if (FALSE && using_cpp_samplers())` to bypass C++
   - Implementing "easy workarounds" that avoid fixing the real issue

### **Current C++ Implementation Status**

- **✅ Working**: Basic MVNormal, MVNormal2, Beta, Weibull, Normal distributions
- **⚠️ Needs Attention**: MVNormal covariance models parameter structure integration
- **🔧 In Progress**: C++ dispatch integration in mixing_distribution_likelihood.R

### **Priority Areas for C++ Development**

1. **High Priority**: Fix MVNormal parameter structure mismatch (see `debug_scripts/mvnormal_cpp_parameter_mismatch_analysis.md`)
2. **Medium Priority**: Complete C++ dispatch for all distributions in mixing_distribution_likelihood.R  
3. **Low Priority**: Optimize C++ implementations for memory usage and speed

### **Current Issue: MVNormal C++ Parameter Structure Mismatch**

**Problem**: Benchmark fails with "values must be length 2, but FUN(X[[1]]) result is length 1" when C++ is enabled
**Root Cause**: `mvnormal_likelihood_wrapper_cpp` function not properly exported and parameter structure handling incorrect
**Analysis**: Complete analysis in `debug_scripts/mvnormal_cpp_parameter_mismatch_analysis.md`

**Required Fixes**:
1. Export `mvnormal_likelihood_wrapper_cpp` function properly
2. Fix multi-cluster parameter handling in `Likelihood.mvnormal` C++ path
3. Complete mvnormal handler in `mixing_distribution_likelihood.R`
4. Test all 9 covariance models with C++ enabled

**Status**: ⚠️ CRITICAL - Blocking benchmark functionality
**Next Steps**: Follow implementation plan in analysis document

## Development Commands

### Package Development (R Console)
- `devtools::test()` - Run the complete test suite
- `devtools::check()` - Full R CMD check (includes tests, documentation, examples)
- `devtools::build()` - Build the package tarball
- `devtools::install()` - Install package locally for testing
- `devtools::document()` - Generate documentation from roxygen2 comments
- `devtools::load_all()` - Load package for interactive development

### Package Development (Bash Terminal)
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::test()"` - Run the complete test suite
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::check()"` - Full R CMD check (includes tests, documentation, examples)  
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::build()"` - Build the package tarball
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::install()"` - Install package locally for testing
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::document()"` - Generate documentation from roxygen2 comments
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::load_all()"` - Load package for interactive development

### Testing
- `testthat::test_check("dirichletprocess")` - Run all tests via testthat (R console)
- `testthat::test_file("tests/testthat/test-filename.R")` - Run specific test file (R console)
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "testthat::test_check('dirichletprocess')"` - Run all tests via testthat (bash terminal)
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "testthat::test_file('tests/testthat/test-filename.R')"` - Run specific test file (bash terminal)

### MVNormal C++ Testing (Known Issues)
- **Issue 1**: `testthat::test_file("tests/testthat/test-mvnormal-cpp-comprehensive.R")` causes R session crashes during cleanup
  - **Root Cause**: testthat/devtools framework has cleanup conflicts with C++ object management
  - **Solution**: Use `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "source('tests/test_mvnormal_cpp_comprehensive_standalone.R')"` for comprehensive C++ testing
- **Issue 2**: MVNormal2 individual C++ functions have parameter format issues [RESOLVED]
  - **Root Cause**: C++ functions expect specific parameter format: `theta <- list(mu_array, sig_array)` not `list(mu=, sig=)`
  - **Error Message**: "Not compatible with requested type: [type=NULL; target=integer]" (misleading error)
  - **Solution**: Use correct theta format for `mvnormal2_likelihood_cpp()`:
    ```r
    # ❌ Wrong: theta <- list(mu = prior_result$mu[1,,1], sig = prior_result$sig[,,1])
    # ✅ Correct: theta <- list(prior_result$mu, prior_result$sig)
    ```
  - **Status**: ✅ RESOLVED - MVNormal2 C++ functions work correctly with proper parameter format
  - **Debug Script**: `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "source('debug_scripts/debug_mvnormal2_type_issues.R')"` demonstrates correct usage

### Benchmarking
- Benchmark scripts are in `benchmark/atime/` directory
- `benchmark-[distribution]-atime.R` - Performance benchmarks for each distribution type
- Results include timing comparisons between R and C++ implementations

### Claude Code R Integration
- **Rscript Path**: `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe"`
- **Usage**: Claude Code can execute R commands via Rscript from bash terminal
- **Examples**:
  - `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(devtools); test()"`
  - `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "source('script.R')"`
- **Benefits**: Enables Claude Code to run R/devtools commands without requiring separate R console interaction

### C++ Development
- Package uses Rcpp and RcppArmadillo for C++ integration
- C++ headers are in `inst/include/`
- C++ source files are in `src/`
- Build system automatically compiles all .cpp files in src/

### Testing Best Practices for C++ Functions

**CRITICAL**: When writing tests for C++ functions, follow these guidelines to ensure environment consistency:

#### C++ Availability Check
```r
# ✅ CORRECT: Check package namespace
cpp_available <- function() {
  if (!using_cpp()) return(FALSE)
  
  pkg_ns <- getNamespace("dirichletprocess")
  required_functions <- c(
    "function_name_cpp",
    "another_function_cpp"
  )
  
  all_exist <- all(sapply(required_functions, function(f) exists(f, where = pkg_ns)))
  return(all_exist)
}

# ❌ WRONG: Check global environment (environment-dependent)
cpp_available <- function() {
  return(using_cpp() && exists("function_name_cpp"))
}
```

#### Accessing C++ Functions in Tests
```r
# ✅ CORRECT: Access functions from package namespace
pkg_ns <- getNamespace("dirichletprocess")
function_name_cpp <- get("function_name_cpp", pkg_ns)
another_function_cpp <- get("another_function_cpp", pkg_ns)

# ❌ WRONG: Direct access (may fail in some environments)
# function_name_cpp() # This may not work consistently
```

#### Test Structure Template
```r
# Enable C++ mode
set_use_cpp(TRUE)

# Check availability with proper namespace checking
cpp_available <- function() {
  if (!using_cpp()) return(FALSE)
  pkg_ns <- getNamespace("dirichletprocess")
  required_functions <- c("your_cpp_functions_here")
  return(all(sapply(required_functions, function(f) exists(f, where = pkg_ns))))
}

# Skip tests if C++ not available
if (!cpp_available()) {
  skip("C++ implementations not available")
}

# Access functions from namespace
pkg_ns <- getNamespace("dirichletprocess")
your_cpp_function <- get("your_cpp_function", pkg_ns)
```

**Why This Matters**:
- Ensures tests work consistently across all development environments
- Follows R package namespace conventions
- Prevents environment-specific test failures
- Maintains proper separation between internal and external functions

## Package Architecture

### Core Object Structure
The package uses S3 classes with a hierarchical inheritance system:
- **`dirichletprocess`**: Main class containing data, clusters, and parameters
- **`mixingDistribution`**: Base class for probability distributions
- **Class hierarchy**: `c("conjugate", "normal", "dirichletprocess")` for method dispatch

### Key Components
- **Data**: `data` field contains the input observations
- **Clusters**: `clusterLabels`, `clusterParameters`, `numberClusters`, `pointsPerCluster`
- **Concentration**: `alpha` parameter with `alphaPriorParameters`
- **Base Distribution**: `mixingDistribution` object defining the component distribution

### MCMC Implementation
The package implements Neal's Algorithm 4 (conjugate) and Algorithm 8 (non-conjugate):
1. **`Fit()`**: Main MCMC sampler using S3 dispatch
2. **`ClusterComponentUpdate()`**: Updates cluster assignments
3. **`ClusterParameterUpdate()`**: Updates cluster parameters
4. **`UpdateAlpha()`**: Updates concentration parameter using West (1992) method

### Distribution Types
- **Conjugate**: Normal-Inverse-Gamma, Normal-Fixed-Variance, Exponential-Gamma, Multivariate Normal
- **Non-Conjugate**: Beta, Weibull, Multivariate Normal 2 (uses Metropolis-Hastings)
- **Hierarchical**: Hierarchical Dirichlet Process variants

### R/C++ Dual Implementation
- **R Implementation**: Pure R code for all algorithms
- **C++ Implementation**: High-performance backend with automatic fallback
- **Interface**: `cpp_interface.R` handles enabling/disabling C++ usage
- **Compatibility**: `can_use_cpp()` checks if C++ implementation is available

Key C++ interface functions:
- `run_mcmc_cpp()`: Main C++ MCMC runner
- `prepare_mixing_dist_params()`: Converts R objects to C++ format
- `enable_cpp_samplers()`: Enables C++ for standard models
- `enable_cpp_hierarchical_samplers()`: Enables C++ for hierarchical models

### Hierarchical vs Standard Models
- **Standard DP**: Single Dirichlet process for clustering
- **Hierarchical DP**: Multiple related DPs sharing cluster parameters
- **Global Parameters**: Shared across groups in hierarchical models
- **Stick-Breaking**: Used in hierarchical implementation

## File Structure

### Core Algorithm Files (R/)
- `dirichletprocess.R`: Main package functions and object creation
- `dirichlet_process_create.R`: Factory functions for DP objects
- `fit.R`: MCMC sampling implementations
- `mixing_distribution.R`: Base distribution class
- `cpp_interface.R`: R/C++ integration layer

### Distribution Implementations (R/)
- `normal_inverse_gamma.R`, `exponential_gamma.R`: Conjugate distributions
- `beta_uniform_gamma.R`, `weibull_uniform_gamma.R`: Non-conjugate distributions
- `mvnormal_normal_wishart.R`: Multivariate normal implementations

### C++ Backend (src/ and inst/include/)
- `mcmc_runner.h/.cpp`: Main MCMC execution engine
- `mixing_distribution_base.h`: Base class for distributions
- `DirichletProcess.h/.cpp`: Core DP implementation
- `*_mixing.h/.cpp`: Distribution-specific implementations

### Testing (tests/testthat/)
- Comprehensive test suite covering all distributions and algorithms
- Tests for R/C++ consistency
- Performance benchmarks in `benchmark/`

## Standard Workflow

1. First think through the problem, read the codebase for relevant files, and write a plan to tasks/todo.md.
2. The plan should have a list of todo items that you can check off as you complete them
3. Before you begin working, check in with me and I will verify the plan.
4. Then, begin working on the todo items, marking them as complete as you go.
5. Every step of the way just give me a high level explanation of what changes you made
6. Finally, add a review section to the todo.md file with a summary of the changes you made and any other relevant information.

## Development Guidelines

### Adding New Distributions
1. Create mixing distribution object with required S3 methods
2. Implement: `Likelihood()`, `PriorDraw()`, `PosteriorDraw()`, `PosteriorParameters()`
3. Add corresponding `DirichletProcess[Distribution]()` constructor
4. Create comprehensive tests
5. Optionally add C++ implementation for performance

### S3 Method Dispatch
The package heavily uses S3 method dispatch based on class inheritance:
- Methods are dispatched on the mixing distribution class
- Conjugate vs non-conjugate algorithms are selected automatically
- Hierarchical models use specialized method implementations

### Performance Considerations
- C++ implementation provides significant speedup for large datasets
- Automatic fallback ensures compatibility if C++ fails
- Vectorized operations are used extensively in R implementation

## Repository Structure

### Documentation and Examples
- **`vignettes/`**: Package documentation with LaTeX vignette and graphics
- **`man/`**: Generated roxygen2 documentation (77 files)
- **`data/`**: Contains `rats.rda` dataset for examples
- **`docs/`**: Generated pkgdown website
- **`papers/`**: Research papers documenting theoretical foundations

### Build Configuration
- **`DESCRIPTION`**: Package metadata, dependencies, and configuration
- **`NAMESPACE`**: Auto-generated package namespace (managed by roxygen2)
- **`_pkgdown.yml`**: Website configuration for documentation
- **`Makevars`/`Makevars.win`**: C++ compilation settings for Unix/Windows

## Important Notes

### Dependencies and Requirements
- **R Version**: R (>= 2.10)
- **C++ Standard**: C++11 required for compilation
- **Key Dependencies**: Rcpp (>= 1.0.11), RcppArmadillo, ggplot2, mvtnorm, gtools
- **System Requirements**: LAPACK, BLAS, optional OpenMP support

### Package URLs
- **Original Repository**: https://github.com/dm13450/dirichletprocess
- **Current Branch**: https://github.com/plato-12/dirichletprocess/tree/cpp-implementation
- **Documentation**: https://dm13450.github.io/dirichletprocess/
- **CRAN**: Available as stable release
- **Issues**: https://github.com/dm13450/dirichletprocess/issues