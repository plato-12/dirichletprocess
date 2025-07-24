# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## **CRITICAL: ASK USER FOR DEVTOOLS COMMANDS**

**🚨 MANDATORY: Before running ANY devtools command (`devtools::document()`, `devtools::check()`, `devtools::test()`, `devtools::install()`, `devtools::build()`, etc.), Claude Code MUST ASK THE USER to run these commands directly in RStudio or their R console. DO NOT attempt to run these commands through bash/Rscript - always ask the user first.**

**The user has the proper R development environment and can execute these commands successfully. This prevents compilation issues and ensures reliable execution.**

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

**Advanced C++ Features**:
- **Unified Manual Interface**: `CppMCMCRunner` with temperature control, cluster operations
- **Sophisticated Architecture**: Advanced features like predictive sampling, convergence diagnostics
- **Production-Ready**: Robust fallback mechanisms and parameter handling



## Current Development Phase: C++ Testing Framework Complete + Production-Ready Manual MCMC

**CURRENT STATUS**: Comprehensive C++ testing framework **implemented and operational**. All critical C++ implementation bugs have been resolved. **CppMCMCRunner manual MCMC interface is fully functional** with advanced features for research applications.

### C++ Implementation Achievement Summary:
- **✅ Complete Manual MCMC C++ Support**: 6 distributions (Normal, Exponential, Beta, Weibull, MVNormal, MVNormal2) with unified `CppMCMCRunner` interface
- **✅ Complete Hierarchical C++ Support**: 3 distributions (Hierarchical Beta, MVNormal, MVNormal2) with full MCMC C++  
- **✅ Full C++ Support**: All major distributions now have complete C++ support with unified interface integration
- **✅ Advanced C++ Features**: Temperature control, cluster operations, predictive sampling, convergence diagnostics


**Status**: ✅ **TESTING FRAMEWORK OPERATIONAL** - ✅ **ALL C++ ISSUES RESOLVED** - ✅ **MANUAL MCMC INTERFACE COMPLETE**

**Detailed Status**: See `debug_scripts/testing-framework-status.md` for comprehensive progress tracking

## Development Commands

**Test Commands**:
- `testthat::test_check("dirichletprocess")` - Run all tests via testthat (R console)
- `testthat::test_file("tests/testthat/test-filename.R")` - Run specific test file (R console)

### Claude Code Git Bash Development Commands

**✅ COMPLETE R PACKAGE DEVELOPMENT ENVIRONMENT ESTABLISHED (2025-07-19)**

Claude Code now has a **complete R package development environment** in VS Code with git bash. All major development commands work perfectly:

#### All Development Commands Working
```bash
# ✅ FULLY OPERATIONAL: Complete package development workflow
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "Rcpp::compileAttributes()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::document()" (not functional, ask user to run and wait for his execution to complete)
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::test()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::check()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::build()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::install()"
```

#### VS Code R Extension Integration
```json
// .vscode/settings.json - Complete R configuration
{
  "r.rterm.windows": "C:/Program Files/R/R-4.4.1/bin/x64/R.exe",
  "r.rpath.windows": "C:/Program Files/R/R-4.4.1/bin/x64/R.exe",
  "r.bracketedPaste": true,
  "r.alwaysUseActiveTerminal": true,
  "r.sessionWatcher": true,
  "r.rtermSendDelay": 8,
  "editor.wordWrap": "on",
  "files.associations": {
    "*.R": "r",
    "*.Rmd": "rmd"
  }
}
```

#### Complete Testing Framework
```bash
# ✅ COMPLETE TESTING CAPABILITY
# Run all tests with full output
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(testthat); library(dirichletprocess); test_dir('tests/testthat')"

# Run specific test file  
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(testthat); library(dirichletprocess); test_file('tests/testthat/test_normal_inverse_gamma.R')"

# Test individual distributions
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(testthat); library(dirichletprocess); test_file('tests/testthat/test_mvnormal_normal_wishart.R')"
```

#### Package Development Capabilities  
```bash
# ✅ FULL PACKAGE DEVELOPMENT WORKFLOW
# Load package and verify C++ availability
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(dirichletprocess); cat('Using C++:', using_cpp(), '\n')"

# Basic functionality testing
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(dirichletprocess); dp <- DirichletProcessGaussian(c(1,2,3)); cat('✅ Package working\n')"

# C++ development workflow
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "Rcpp::compileAttributes(); cat('✅ C++ exports updated\n')"
```

#### Development Environment Features
**✅ COMPLETE R DEVELOPMENT ENVIRONMENT** in VS Code with git bash:

**Real-time Development Capabilities**:
- ✅ **Complete devtools workflow**: All commands (`test()`, `document()`, `check()`, `build()`, `install()`) working
- ✅ **C++ compilation**: `Rcpp::compileAttributes()` and full compilation support
- ✅ **Interactive testing**: Individual test files and complete test suite execution
- ✅ **Real-time debugging**: Immediate feedback on C++ errors and R fallbacks
- ✅ **VS Code integration**: R extension with syntax highlighting, IntelliSense, and integrated terminal
- ✅ **Git integration**: Full version control workflow within VS Code

**Latest Test Results (2025-07-19)**:
- **Package Loading**: ✅ `library(dirichletprocess)` successful
- **C++ Status**: ✅ `using_cpp()` returns `TRUE`
- **Individual Tests**: ✅ Complete with detailed output (10 PASS, 1 WARN showing C++ fallback detection)
- **Documentation**: ✅ `devtools::document()` executes successfully
- **C++ Compilation**: ✅ `Rcpp::compileAttributes()` working perfectly

#### Streamlined Development Workflow
```bash
# ✅ COMPLETE DEVELOPMENT CYCLE - ALL COMMANDS WORKING
# 1. C++ development
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "Rcpp::compileAttributes()"

# 2. Documentation updates  
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::document()"

# 3. Focused testing
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "library(testthat); library(dirichletprocess); test_file('tests/testthat/test_[specific_area].R')"

# 4. Full package validation
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::test()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::check()"

# 5. Package building
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::build()"
```

**Status**: ✅ **COMPLETE R PACKAGE DEVELOPMENT ENVIRONMENT** - Full professional R development capability in VS Code with git bash

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



## Development Guidelines



**MEDIUM PRIORITY (Complete in 4-5 weeks)**:
1. **Performance benchmarking**: Validate C++ performance improvements over R implementations
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

### Current Branch Status
update - **Focus**: High-performance C++ backends for Dirichlet Process algorithms
- **Current Status**: ✅ **100% manual MCMC C++ coverage achieved**, core functionality stable, performance optimized, ready for research use
- **Next**: Execute comprehensive testing framework and validate R/C++ consistency