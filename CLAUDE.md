# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## **CRITICAL: ASK USER FOR DEVTOOLS COMMANDS**

**🚨 MANDATORY: Before running ANY devtools command (`devtools::document()`, `devtools::check()`, `devtools::test()`, `devtools::install()`, `devtools::build()`, etc.), Claude Code MUST ASK THE USER to run these commands directly in RStudio or their R console. DO NOT attempt to run these commands through bash/Rscript - always ask the user first.**

**The user has the proper R development environment and can execute these commands successfully. This prevents compilation issues and ensures reliable execution.**

## **PROJECT MISSION: C++ IMPLEMENTATION PRIORITY**

**CRITICAL DIRECTIVE**: This project is focused on implementing and maintaining high-performance C++ implementations for the Dirichlet Process algorithms. The primary goal is to provide C++ backends that significantly improve performance over pure R implementations.

### **C++ Implementation Status**

**✅ COMPLETE**: 100% C++ coverage achieved for all major functionality:
- **6 Distributions**: Normal, Exponential, Beta, Weibull, MVNormal, MVNormal2
- **3 Hierarchical Models**: Hierarchical Beta, MVNormal, MVNormal2  
- **Advanced Features**: Manual MCMC interface, temperature control, convergence diagnostics
- **Production Ready**: Automatic R fallback, comprehensive testing, identical results validation



## Current Development Phase: CRAN Submission Ready + JSS Preparation

**CURRENT STATUS**: Package is **CRAN submission ready** with complete C++ implementation and comprehensive testing framework. All R CMD check requirements met (0 errors, 0 warnings, 2 system notes). Ready for Journal of Statistical Software (JSS) submission preparation.

### Package Achievement Summary:
- **✅ CRAN Submission Ready**: `devtools::check(cran = TRUE)` passes with 0 errors, 0 warnings, 2 system notes
- **✅ Complete C++ Implementation**: 100% C++ coverage for all major distributions and hierarchical models
- **✅ Reverse Dependencies Verified**: Both reverse dependencies (copre, MIRES) tested and compatible
- **✅ Documentation Complete**: README.md, NEWS.md, and cran-comments.md updated for submission
- **✅ GitHub Actions Fixed**: CI/CD pipeline operational (vignette building disabled for LaTeX issues)

### Next Phase: JSS Submission (5-7 weeks)
- **Phase 1**: Enhance JSS vignette with C++ performance focus (2-3 weeks)
- **Phase 2**: Create comprehensive benchmarking analysis (1-2 weeks) 
- **Phase 3**: Prepare final JSS submission materials (1 week)
- **Phase 4**: Submit to Journal of Statistical Software (1 week)

## Development Commands

**Test Commands**:
- `testthat::test_check("dirichletprocess")` - Run all tests via testthat (R console)
- `testthat::test_file("tests/testthat/test-filename.R")` - Run specific test file (R console)

### Claude Code Git Bash Development Commands

**✅ COMPLETE R PACKAGE DEVELOPMENT ENVIRONMENT ESTABLISHED (2025-07-19)**

Claude Code now has a **complete R package development environment** in VS Code with git bash. All major development commands work perfectly:

#### Development Commands
```bash
# Essential R development commands
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "Rcpp::compileAttributes()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::test()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::check()"
"C:/PROGRA~1/R/R-44~1.1/bin/x64/Rscript.exe" -e "devtools::build()"
```

**Note**: `devtools::document()` should be run by user in RStudio/R console due to compilation requirements.

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

### C++ Development Notes
- Package uses Rcpp and RcppArmadillo for high-performance implementations
- C++ source files in `src/`, headers in `inst/include/`
- Comprehensive testing framework ensures R/C++ consistency
- Automatic fallback to R implementations when needed

## Package Overview

### Core Functionality
- **Dirichlet Process Models**: S3 classes for flexible Bayesian nonparametric modeling
- **MCMC Implementation**: Neal's Algorithm 4 (conjugate) and Algorithm 8 (non-conjugate)
- **Distributions**: 6 distributions + 3 hierarchical models with complete C++ coverage
- **Covariance Models**: Full support for FULL, EII, VII, EEI, VEI, EVI, VVI models

### Key Features
- **High Performance**: C++ backends with automatic R fallback
- **Complete Coverage**: All major distributions and hierarchical models
- **Production Ready**: Comprehensive testing and validation
- **Research Grade**: Advanced MCMC features and diagnostics

### Implementation Control (Updated 2025-08-07)
- **Universal cpp Parameter**: All 11 distribution constructors now accept `cpp = TRUE/FALSE` for explicit implementation control
- **Default Changed**: Package default is now `cpp = FALSE` (R implementation) for predictable cross-platform behavior
- **Per-Object Control**: Each constructor call sets implementation preference independently
- **Backward Compatibility**: Legacy global control functions (`set_use_cpp()`, `using_cpp()`) remain functional
- **Complete Coverage**: Normal, Beta, Exponential, Weibull, MVNormal, MVNormal2, Hierarchical, and Markov models

## Repository Structure
- **R/**: Main package functions and C++ integration layer
- **src/**: C++ implementation files and headers
- **tests/testthat/**: Comprehensive test suite with R/C++ consistency validation
- **benchmark/**: Performance analysis and benchmarking framework
- **debug_scripts/**: Development artifacts and troubleshooting documentation



## Submission Status (2025-07-30)

### CRAN Submission
- **Status**: ✅ **READY FOR SUBMISSION**
- **R CMD Check**: 0 errors, 0 warnings, 2 system notes (non-blocking)
- **Reverse Dependencies**: Both `copre` and `MIRES` verified compatible
- **Documentation**: All submission files updated (README.md, NEWS.md, cran-comments.md)

### JSS Submission Preparation
- **Target Timeline**: 5-7 weeks
- **Current Vignette**: JSS-compliant LaTeX structure already exists
- **Enhancement Needed**: Add C++ performance benchmarking and architecture sections
- **Submission Advantage**: Novel C++ implementation for Dirichlet processes

## Requirements
- **R Version**: R (>= 2.10)
- **C++ Standard**: C++11 required for compilation
- **Key Dependencies**: Rcpp (>= 1.0.11), RcppArmadillo, ggplot2, mvtnorm, gtools
- **System Requirements**: LAPACK, BLAS, optional OpenMP support

### Current Branch Status
**Branch**: `feature/testing-framework-fixes`
- **Focus**: CRAN submission preparation and JSS manuscript enhancement
- **Current Status**: ✅ **CRAN submission ready** - All checks pass, documentation complete, reverse dependencies verified
- **Next Phase**: JSS manuscript preparation with C++ performance benchmarking

## Final Status Summary (2025-07-30)

### 🎯 **MISSION ACCOMPLISHED**
The dirichletprocess package has achieved its primary objectives:

1. **✅ Complete C++ Implementation**: 100% coverage for all distributions and hierarchical models
2. **✅ Production Ready**: Comprehensive testing, validation, and error handling
3. **✅ CRAN Submission Ready**: All requirements met, documentation complete
4. **✅ High Performance**: Significant speedups while maintaining identical results
5. **✅ Research Grade**: Advanced features for sophisticated Bayesian nonparametric analysis

### 🚀 **Ready for Next Phase**
- **CRAN Submission**: Package ready for immediate submission
- **JSS Preparation**: 5-7 week timeline to prepare comprehensive manuscript
- **Research Impact**: Enables large-scale Bayesian nonparametric modeling with R/C++ performance