# MVNormal2 CppMCMCRunner Architecture Analysis

## Core Issues Identified

### 1. **Algorithm Mismatch in Testing Framework**
The main architectural problem is in `helper-testing.R` lines 176-181:
```r
# C++ implementation
set_use_cpp(TRUE)
# Enable C++ samplers
enable_cpp_samplers(TRUE)
enable_cpp_hierarchical_samplers(TRUE)
dp_cpp <- create_dp_object(distribution_type, test_data)
dp_cpp <- Fit(dp_cpp, its = iterations)
```

**Problem**: The testing framework calls `enable_cpp_samplers(TRUE)` and `enable_cpp_hierarchical_samplers(TRUE)` but these functions don't exist with parameters. They are just status check functions.

**Root Cause**: The helper functions assume these are setter functions, but they are actually status checkers:
- `enable_cpp_samplers()` returns `TRUE/FALSE` (line 260 in cpp_interface.R)
- `enable_cpp_hierarchical_samplers()` doesn't exist

### 2. **Parameter Format Inconsistency** 
MVNormal2 expects parameters in different formats between R and C++:
- **R implementation**: Uses named list with `mu` and `sig` arrays
- **C++ implementation**: Uses flattened parameter vectors (d elements for mu + d*(d+1)/2 for Sigma upper triangle)

### 3. **CppMCMCRunner vs Direct Function Call Confusion**
The current architecture has two different C++ interfaces:
1. **Unified CppMCMCRunner**: `run_mcmc_cpp()` through MCMCRunner class
2. **Individual function calls**: Direct calls to `mvnormal2_likelihood_cpp()`, etc.

**Problem**: The MVNormal2 R code tries to use individual function calls (line 21 in mvnormal_semi_conjugate.R), but the consistency tests try to use the unified interface.

### 4. **Missing Enable Functions**
The testing framework expects enable functions that don't exist:
```r
enable_cpp_samplers(TRUE)           # This function doesn't take parameters
enable_cpp_hierarchical_samplers(TRUE)  # This function doesn't exist at all
```

### 5. **Timeout Due to Inconsistent Random Number Generation**
The consistency tests run 5 runs of 50-200 iterations each, but there may be RNG seed issues between R and C++ implementations causing very slow convergence or infinite loops in some cases.

## Architectural Solutions Required

### 1. **Fix Testing Framework Enable Functions**
Create proper enable/disable functions in `cpp_interface.R`:
```r
enable_cpp_samplers <- function(enable = TRUE) {
  # Actual implementation to enable/disable C++ samplers
  options(dirichletprocess.force_cpp_samplers = enable)
}

enable_cpp_hierarchical_samplers <- function(enable = TRUE) {
  options(dirichletprocess.force_cpp_hierarchical = enable)
}
```

### 2. **Standardize Parameter Format**
Ensure MVNormal2 uses consistent parameter format between R and C++ by:
- Converting R format (mu/sig arrays) to C++ format (flattened vectors) in interface
- Or modifying C++ to accept R format directly

### 3. **Unify CppMCMCRunner Integration**
Choose one approach:
- **Option A**: Use unified MCMCRunner for all distributions including MVNormal2
- **Option B**: Use direct function calls for all distributions
- **Option C**: Hybrid approach with clear documentation

### 4. **Fix Algorithm Selection**
Ensure MVNormal2 (non-conjugate) uses Algorithm 8 in both R and C++ implementations consistently.

## Recommended Fix Priority

1. **HIGH**: Fix missing enable functions (causes immediate test failures)
2. **HIGH**: Fix parameter format consistency (causes result differences)  
3. **MEDIUM**: Unify CppMCMCRunner architecture
4. **LOW**: Optimize performance and timeouts

## Current Status
- Basic MVNormal2 C++ functionality works (confirmed with simple test)
- Consistency testing framework has architectural mismatches
- Need to choose unified approach for C++ integration