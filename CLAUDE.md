# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Package Development (R Console)
- `devtools::test()` - Run the complete test suite
- `devtools::check()` - Full R CMD check (includes tests, documentation, examples)
- `devtools::build()` - Build the package tarball
- `devtools::install()` - Install package locally for testing
- `devtools::document()` - Generate documentation from roxygen2 comments
- `devtools::load_all()` - Load package for interactive development

### Package Development (Bash Terminal)
- `R -e "devtools::test()"` - Run the complete test suite
- `R -e "devtools::check()"` - Full R CMD check (includes tests, documentation, examples)
- `R -e "devtools::build()"` - Build the package tarball
- `R -e "devtools::install()"` - Install package locally for testing
- `R -e "devtools::document()"` - Generate documentation from roxygen2 comments
- `R -e "devtools::load_all()"` - Load package for interactive development

### Testing
- `testthat::test_check("dirichletprocess")` - Run all tests via testthat (R console)
- `testthat::test_file("tests/testthat/test-filename.R")` - Run specific test file (R console)
- `R -e "testthat::test_check('dirichletprocess')"` - Run all tests via testthat (bash terminal)
- `R -e "testthat::test_file('tests/testthat/test-filename.R')"` - Run specific test file (bash terminal)

### C++ Development
- Package uses Rcpp and RcppArmadillo for C++ integration
- C++ headers are in `inst/include/`
- C++ source files are in `src/`
- Build system automatically compiles all .cpp files in src/

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
6. Folder papers/ contain the research papers relevant to the algorithms used in the package.

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