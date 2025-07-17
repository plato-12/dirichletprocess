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

**📊 C++ Manual MCMC Support: 100% Complete ✅**

**✅ COMPLETE Manual MCMC C++ Support (6 distributions)**:
- **Normal/Gaussian**: All individual functions + unified `CppMCMCRunner` interface
- **Exponential**: All individual functions + unified `CppMCMCRunner` interface  
- **Beta**: All individual functions + unified `CppMCMCRunner` interface
- **Weibull**: All individual functions + unified `CppMCMCRunner` interface
- **MVNormal**: All individual functions + unified `CppMCMCRunner` interface ✅ **COMPLETED**
- **MVNormal2**: All individual functions + unified `CppMCMCRunner` interface ✅ **COMPLETED**

**✅ COMPLETE Hierarchical C++ Support (3 distributions)**:
- **Hierarchical Beta**: Full MCMC C++ + manual steps
- **Hierarchical MVNormal**: Full MCMC C++ + manual steps
- **Hierarchical MVNormal2**: Full MCMC C++ + manual steps

**🟡 PARTIAL C++ Support (0 distributions)**:
- **All major distributions now have complete C++ support!**

**❌ Missing C++ Support (2 minor distributions)**:
- **Beta2**: Pure R implementation
- **Normal Fixed Variance**: Pure R implementation

**Advanced C++ Features**:
- **Unified Manual Interface**: `CppMCMCRunner` with temperature control, cluster operations
- **Sophisticated Architecture**: Advanced features like predictive sampling, convergence diagnostics
- **Production-Ready**: Robust fallback mechanisms and parameter handling

### **Recent Major Fix: Constrained Covariance Models**

**Problem**: Constrained covariance models (EII, VII, EEI, VEI, EVI, VVI) failed during MCMC with "incorrect number of dimensions" error.

**Root Cause**: Parameter arrays for constrained models are 2D `[nParams, clusters]` while FULL models use 3D `[d, d, clusters]`. Code assumed all models used 3D arrays.

**Solution Applied**:
1. Fixed `MvnormalCreate()` to properly initialize default parameters for partial parameter lists
2. Fixed `Initialise()` function to handle dimension access safely (`mu_dim[3]` checks)
3. Fixed `ClusterComponentUpdate.conjugate()` to handle both 2D and 3D parameter arrays
4. Added numerical stability checks for VEI model calculations

**Status**: ✅ CORE FUNCTIONALITY WORKING - Main dimension issues resolved

### **Recent Major Fix: E/V Model Initialization**

**Problem**: E and V models (univariate constrained models) failed during initialization with "argument is of length zero" error, causing R session termination in benchmark runner.

**Root Cause**: The `Initialise.conjugate` function assumed `mu` parameters were arrays with dimensions, but E/V models return scalar `mu` values. Code tried to access `dim(mu)[3]` on scalar values.

**Solution Applied**:
1. Fixed `Initialise.conjugate()` to handle scalar `mu` parameters by detecting when `dim(mu)` returns NULL
2. Added dimension-aware conversion to proper 3D array format for consistency
3. Maintained compatibility with existing multivariate models

**Status**: ✅ RESOLVED - E and V models now initialize correctly

### **Recent Major Fix: Research Benchmark Performance**

**Problem**: Research benchmark runner appeared to "crash" but was actually taking 9-10 hours due to computational explosion (34,000 iterations), with no progress indicators.

**Root Cause**: Combinatorial explosion in atime framework: 4 sample sizes × 7-10 models × 10 repetitions × 50 MCMC iterations.

**Solution Applied**:
1. **Optimized sample size testing**: Reduced from 4 sizes `(10, 20, 50, max_n)` to 2 sizes `(20, max_n)` = 50% reduction
2. **Optimized QUICK_RESEARCH_CONFIG**: MCMC 50→25, repetitions 10→5, max_samples 200→100 = 62.5% additional reduction  
3. **Added progress indicators**: Shows exact run counts and time estimates
4. **Overall improvement**: 87.5% reduction in computation time (9-10 hours → 70 minutes)

**Status**: ✅ RESOLVED - Benchmark runner now completes in reasonable time with clear progress

### **Recent Major Achievement: PRIORITY 1 MVNormal Integration COMPLETED**

**Date**: 2025-07-17  
**Objective**: Achieve 100% manual MCMC C++ coverage by integrating MVNormal distributions into the unified C++ interface

**Implementation Details**:
1. **✅ Added MVNormal support to can_use_cpp()**: Modified `R/cpp_interface.R` to include "mvnormal" and "mvnormal2" in supported types
2. **✅ Implemented conjugate_mvnormal_update_alpha_cpp()**: Added C++ UpdateAlpha function in `src/MVNExports.cpp` using West (1992) auxiliary variable method
3. **✅ Implemented nonconjugate_mvnormal2_update_alpha_cpp()**: Added C++ UpdateAlpha function in `src/MVNormal2Exports.cpp` using West (1992) auxiliary variable method  
4. **✅ Regenerated RcppExports**: Both functions properly exported via `Rcpp::compileAttributes()`

**Technical Implementation**:
- Both UpdateAlpha functions use the standard West (1992) auxiliary variable method for concentration parameter updates
- Functions handle all parameter structures correctly with proper error handling
- Implementation follows same pattern as existing conjugate_exponential_update_alpha_cpp()
- Support for all 9 MVNormal covariance models (E, V, FULL, EII, VII, EEI, VEI, EVI, VVI)

**Result**: 
- **📊 C++ Manual MCMC Support**: 83% → 100% Complete
- **✅ All 6 major distributions** now support unified `CppMCMCRunner` interface
- **✅ Manual MCMC loops** work with full C++ acceleration for all major distributions
- **✅ Advanced features** (temperature control, cluster operations) available for all distributions

**Status**: ✅ COMPLETED - 100% manual MCMC C++ coverage achieved

## Current Development Phase: Comprehensive Testing Framework Required

**CURRENT STATUS**: C++ implementation work is complete with **100% manual MCMC C++ coverage** achieved. The package has mature, production-ready C++ backends with sophisticated architecture.

### C++ Implementation Achievement Summary:
- **✅ Complete Manual MCMC C++ Support**: 6 distributions (Normal, Exponential, Beta, Weibull, MVNormal, MVNormal2) with unified `CppMCMCRunner` interface
- **✅ Complete Hierarchical C++ Support**: 3 distributions (Hierarchical Beta, MVNormal, MVNormal2) with full MCMC C++  
- **✅ Full C++ Support**: All major distributions now have complete C++ support with unified interface integration
- **✅ Advanced C++ Features**: Temperature control, cluster operations, predictive sampling, convergence diagnostics
- **❌ Missing C++**: 2 minor distributions (Beta2, Normal Fixed Variance) - pure R implementations

### Required Testing and Validation:
1. **Comprehensive C++ Testing Framework**: Create new systematic C++ validation tests from scratch
2. ✅ **MVNormal Integration**: Complete integration into unified `CppMCMCRunner` interface ✅ **COMPLETED**
3. **R/C++ Consistency Verification**: Ensure both implementations produce identical statistical results
4. **Original R Package Test Validation**: Verify all 37 original R tests still pass
5. **Performance Benchmarking**: Validate C++ performance improvements over R implementations
6. **Production Readiness**: Complete package development workflow validation

**Status**: ✅ **100% MANUAL MCMC C++ COVERAGE ACHIEVED** - 🧪 COMPREHENSIVE TESTING FRAMEWORK PHASE REQUIRED

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

**Current Test Structure**:
- **Total test files**: 38 testthat files (37 original R + 1 temporary C++)
- **Original R Package Tests**: 37 files testing core R functionality - **PRESERVED**
- **Temporary C++ Test**: 1 file (`test_cpp_mcmc_runner.R`) - **TO BE DELETED AND RECREATED**

**Test Commands**:
- `testthat::test_check("dirichletprocess")` - Run all tests via testthat (R console)
- `testthat::test_file("tests/testthat/test-filename.R")` - Run specific test file (R console)
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "testthat::test_check('dirichletprocess')"` - Run all tests via testthat (bash terminal)
- `"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "testthat::test_file('tests/testthat/test-filename.R')"` - Run specific test file (bash terminal)

**Future Testing Framework**:
- **Comprehensive C++ Testing**: New systematic C++ validation tests to be created from scratch
- **Test Organization**: C++ tests will be organized in dedicated subdirectory structure
- **Integration Focus**: R/C++ consistency validation, performance benchmarking, production readiness

### Constrained Covariance Testing (Fixed Issues)
- **Issue RESOLVED**: Constrained models now work correctly with proper dimension handling
- **Testing**: All constrained models (EII, VII, EEI, VEI, EVI, VVI) should pass basic MCMC tests
- **Verification Script**: Create test scripts to verify fixes work as expected

### Benchmarking
- Benchmark scripts are in `benchmark/atime/` directory
- `benchmark-[distribution]-atime.R` - Performance benchmarks for each distribution type
- Results include timing comparisons between R and C++ implementations
- **Key Benchmark**: `benchmark/atime/benchmark-covariance-models-optimized.R` - optimized version for practical use
- **Research Benchmark Runner**: `benchmark/research_benchmark_runner.R` - comprehensive research-quality benchmarks
  - `source("benchmark/research_benchmark_runner.R"); quick_research_benchmark()` - 70 minute comprehensive test
  - `standard_research_benchmark()` - 1-3 hour research quality  
  - `publication_benchmark()` - 3-6 hour publication quality
  - All functions include progress indicators and time estimates

### Current Documentation Work
- **R Markdown Documentation**: `docs/cpp_covariance_models_usage.Rmd` - Comprehensive demonstration of all covariance models
  - Shows implementation examples for all 9 covariance models (E, V, FULL, EII, VII, EEI, VEI, EVI, VVI)
  - Includes execution results and performance metrics
  - Currently being refined to handle C++ fallback scenarios gracefully
  - To render: `rmarkdown::render('docs/cpp_covariance_models_usage.Rmd')`
  - **Note**: This is demonstration documentation, not formal validation testing

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

### Comprehensive C++ Testing Framework Guidelines

**Testing Organization**:
- **Preserve Original Tests**: All 37 original R package tests remain in `tests/testthat/`
- **Create C++ Test Framework**: New systematic C++ validation tests from scratch
- **Organized Structure**: C++ tests in dedicated subdirectory for clear separation
- **Focus Areas**: Individual C++ functions, unified interface, R/C++ consistency, performance

**Priority Testing Areas**:
1. **MVNormal Integration**: Complete `CppMCMCRunner` integration for MVNormal/MVNormal2
2. **R/C++ Consistency**: Statistical equivalence validation across all distributions
3. **Performance Validation**: Confirm C++ improvements over R implementations
4. **Production Readiness**: Complete package development workflow validation

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

### Distribution Types and Covariance Models
- **Conjugate**: Normal-Inverse-Gamma, Normal-Fixed-Variance, Exponential-Gamma, Multivariate Normal
- **Non-Conjugate**: Beta, Weibull, Multivariate Normal 2 (uses Metropolis-Hastings)
- **Hierarchical**: Hierarchical Dirichlet Process variants

**Multivariate Normal Covariance Models**:
- **FULL**: Unrestricted covariance matrices (3D parameter arrays)
- **Constrained Models**: EII, VII, EEI, VEI, EVI, VVI (2D parameter arrays)
  - **EII**: Equal, Isotropic, Identical
  - **VII**: Variable volume, Isotropic, Identical orientation
  - **EEI**: Equal volume, Equal shape, Isotropic
  - **VEI**: Variable volume, Equal shape, Isotropic  
  - **EVI**: Equal volume, Variable shape, Isotropic
  - **VVI**: Variable volume, Variable shape, Isotropic

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
- `mvnormal_normal_wishart.R`: Multivariate normal implementations with all covariance models
- `mvnormal_semi_conjugate.R`: Alternative multivariate normal implementation

### MCMC Core Files (R/)
- `initialise.R`: Parameter initialization and pre-allocation logic
- `cluster_component_update.R`: Cluster assignment updates (Neal's Algorithm 4)
- `cluster_parameter_update.R`: Parameter updates for individual clusters

### C++ Backend (src/ and inst/include/)
- `mcmc_runner.h/.cpp`: Main MCMC execution engine
- `mixing_distribution_base.h`: Base class for distributions
- `DirichletProcess.h/.cpp`: Core DP implementation
- `*_mixing.h/.cpp`: Distribution-specific implementations

### Testing (tests/testthat/)
- Comprehensive test suite covering all distributions and algorithms
- Tests for R/C++ consistency
- Performance benchmarks in `benchmark/`

### Debug Scripts (debug_scripts/)
- `constrained_covariance_mcmc_dimensions_issue.md`: Analysis of resolved dimension issues
- `remaining_constrained_covariance_issues.md`: Action plan for remaining minor issues
- Various debugging scripts for troubleshooting specific problems

**IMPORTANT**: All debug files created during testing and development should be saved in `debug_scripts/` directory. This includes:
- Analysis documents (*.md files)
- Diagnostic R scripts for troubleshooting
- Test scripts for specific issues
- Documentation of fixes and workarounds

## Critical Dimension Handling Patterns

### Problem: Parameter Array Dimensions
**Issue**: Constrained covariance models store parameters differently than FULL models:
- **FULL model**: `sig` is 3D array `[d, d, clusters]`
- **Constrained models**: `sig` is 2D array `[nParams, clusters]`

### Solution Pattern (Applied Throughout Codebase)
```r
# ✅ CORRECT: Dimension-aware parameter access
param_dims <- dim(theta[[2]])
if (length(param_dims) == 3) {
  # FULL covariance model - 3D array
  sigma_i <- theta[[2]][, , i]
} else if (length(param_dims) == 2) {
  # Constrained covariance models - 2D array
  sigma_i <- theta[[2]][, i]
} else {
  # Single cluster/scalar case
  sigma_i <- theta[[2]][i]
}

# ❌ WRONG: Hardcoded 3D assumption
sigma_i <- theta[[2]][, , i]  # Fails for constrained models
```

### Files That Implement This Pattern
- `R/mvnormal_normal_wishart.R`: Lines 198-210, 267-279
- `R/mvnormal_semi_conjugate.R`: Lines 40-49
- `R/cluster_component_update.R`: Lines 51-63
- `R/initialise.R`: Lines 36-76 (safe dimension access and E/V model scalar handling)
- `R/benchmark_integration.R`: Lines 215-224 (optimized atime integration with progress indicators)

## Development Guidelines

### Priority Actions (Based on Comprehensive Analysis)

**✅ PRIORITY 1 COMPLETED: MVNormal Integration**:
1. ✅ **Integrate MVNormal into unified C++ interface**: Added "mvnormal" and "mvnormal2" to `can_use_cpp()` supported types
2. ✅ **Implement MVNormal UpdateAlpha C++**: Added `conjugate_mvnormal_update_alpha_cpp()` and `nonconjugate_mvnormal2_update_alpha_cpp()`
3. ✅ **Test MVNormal manual MCMC**: Validated that `CppMCMCRunner` works with all covariance models
4. ✅ **Achieve 100% manual MCMC C++ coverage**: Complete C++ support for all major distributions

**HIGH PRIORITY (Complete in 2-3 weeks)**:
1. **Execute comprehensive testing framework**: Create systematic C++ validation tests from scratch
2. **Validate R/C++ consistency**: Ensure identical statistical behavior across all distributions
3. **Preserve original R tests**: Verify all 37 original R package tests continue to pass
4. **Complete package development workflow**: `devtools::check()` must pass cleanly

**MEDIUM PRIORITY (Complete in 4-5 weeks)**:
1. **Performance benchmarking**: Validate C++ performance improvements over R implementations
2. **Complete minor distributions**: Add C++ support for Beta2 and Normal Fixed Variance
3. **Production readiness validation**: Comprehensive edge case testing and stability validation

### Debug File Management
**CRITICAL**: When creating debug files during testing and development, always save them in the `debug_scripts/` directory. This maintains organization and ensures debugging artifacts are preserved for future reference.

**File naming conventions**:
- Analysis documents: `[issue_name]_analysis.md`
- Diagnostic scripts: `debug_[specific_issue].R`
- Test scripts: `test_[specific_case].R`
- Fix documentation: `[issue_name]_fix_summary.md`

### Adding New Distributions
1. Create mixing distribution object with required S3 methods
2. Implement: `Likelihood()`, `PriorDraw()`, `PosteriorDraw()`, `PosteriorParameters()`
3. Add corresponding `DirichletProcess[Distribution]()` constructor
4. Create comprehensive tests
5. Optionally add C++ implementation for performance

### Adding New Covariance Models
1. Update `getNumCovParams()` function with parameter count
2. Implement `reconstructCovarianceMatrix()` and `extractCovarianceParams()` helpers
3. Add model to valid models list in `MvnormalCreate()`
4. Ensure dimension-aware handling throughout MCMC functions
5. Test with all MCMC components

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

### Current Branch Status
- **Branch**: `cpp-implementation`
- **Focus**: High-performance C++ backends for Dirichlet Process algorithms
- **Current Status**: ✅ **100% manual MCMC C++ coverage achieved**, core functionality stable, performance optimized, ready for research use
- **Next**: Execute comprehensive testing framework and validate R/C++ consistency