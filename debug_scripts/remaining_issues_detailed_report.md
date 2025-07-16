# Remaining Issues: Detailed Report and Action Plan

## Current Status Overview

✅ **RESOLVED**: Core dimension access issues in R code  
⚠️ **REMAINING**: Two separate issues preventing full constrained covariance functionality

---

## Issue 1: C++ Implementation Matrix Dimension Bugs

### **Problem Description**

**Error**: `"addition: incompatible matrix dimensions: 0x1 and 2x1"`  
**Context**: Occurs when C++ samplers are enabled (`enable_cpp_samplers()`)  
**Location**: C++ backend implementation for MVNormal constrained models  
**Trigger**: During initialization or MCMC operations with constrained covariance models

### **Evidence**

```r
# With C++ samplers enabled
library(dirichletprocess)
enable_cpp_samplers()
set.seed(123)
x <- matrix(rnorm(20), ncol = 2)
dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = "EII")))
dp <- Initialise(dp)
# Error: addition: incompatible matrix dimensions: 0x1 and 2x1
```

**Working**: FULL model likely works with C++ (needs verification)  
**Failing**: All constrained models (EII, VII, EEI, VEI, EVI, VVI)

### **Root Cause Analysis**

1. **C++ Parameter Handling**: The C++ implementation may not correctly handle the 2D parameter arrays used by constrained models
2. **Matrix Operations**: C++ matrix operations are expecting different dimensions than what constrained models provide
3. **Parameter Conversion**: The interface between R and C++ may not correctly convert constrained model parameters

### **Likely Locations**

Based on available C++ functions:
- `mvnormal_posterior_draw_cpp()` - Primary suspect for parameter sampling
- `mvnormal_posterior_parameters_cpp()` - May have dimension assumptions
- `conjugate_mvnormal_cluster_parameter_update_cpp()` - Cluster parameter handling
- `conjugate_mvnormal_cluster_component_update_cpp()` - Component updates

### **Investigation Steps**

#### **Phase 1: Isolate the Problem**
1. **Test C++ vs R Dispatch**:
   ```r
   # Test if FULL model works with C++
   enable_cpp_samplers()
   test_full_model_cpp()
   
   # Test if disabling C++ resolves constrained model issues
   set_use_cpp(FALSE)
   test_constrained_models_r_only()
   ```

2. **Identify Exact Failure Point**:
   ```r
   # Add debug prints to isolate where the error occurs
   options(error = function() {
     cat("Call stack at error:\n")
     traceback()
   })
   ```

3. **Test Individual C++ Functions**:
   ```r
   # Test each C++ function individually
   test_mvnormal_posterior_draw_cpp()
   test_mvnormal_posterior_parameters_cpp()
   ```

#### **Phase 2: Analyze C++ Implementation**
1. **Examine C++ Source Code**:
   ```bash
   # Look for matrix dimension assumptions in C++ code
   grep -r "matrix.*dim" src/
   grep -r "0x1\|2x1" src/
   ```

2. **Check Parameter Conversion**:
   ```r
   # Examine how R parameters are converted to C++ format
   analyze_parameter_conversion()
   ```

3. **Compare FULL vs Constrained Handling**:
   ```r
   # Compare parameter structures between models
   compare_parameter_structures()
   ```

#### **Phase 3: Fix Implementation**
1. **Update C++ Parameter Handling**:
   - Modify C++ functions to handle 2D parameter arrays
   - Update matrix operations to work with constrained model dimensions
   - Fix parameter conversion between R and C++

2. **Test Fixes**:
   ```r
   # Comprehensive testing of all constrained models
   test_all_constrained_models_cpp()
   ```

### **Implementation Plan**

#### **Step 1: Diagnostic Analysis**
- [ ] Create test script to isolate C++ vs R behavior
- [ ] Identify exact failure point in C++ call stack
- [ ] Analyze parameter structures for FULL vs constrained models
- [ ] Document exact matrix dimension mismatches

#### **Step 2: C++ Code Investigation**
- [ ] Examine C++ source files for matrix dimension assumptions
- [ ] Identify parameter conversion logic between R and C++
- [ ] Map out C++ function call flow for constrained models
- [ ] Compare with working FULL model implementation

#### **Step 3: C++ Implementation Fix**
- [ ] Update C++ functions to handle 2D parameter arrays
- [ ] Fix matrix operations for constrained model dimensions
- [ ] Update parameter conversion logic
- [ ] Test individual C++ functions with constrained models

#### **Step 4: Integration Testing**
- [ ] Test all constrained models with C++ enabled
- [ ] Verify no regression in FULL model functionality
- [ ] Performance testing of C++ vs R implementations

---

## Issue 2: R Fallback Numerical Stability Issues

### **Problem Description**

**Error**: `"BLAS/LAPACK routine 'DPOTRFDSTEBZDSPEVXDS' gave error code -4"`  
**Context**: Occurs when C++ is disabled and R fallback is used  
**Location**: `PosteriorDraw.mvnormal()` function, specifically in `rWishart()` call  
**Trigger**: During parameter sampling for any MVNormal model (including FULL)

### **Evidence**

```r
# With C++ disabled
library(dirichletprocess)
set_use_cpp(FALSE)
set.seed(123)
x <- matrix(rnorm(20), ncol = 2)
dp <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = "FULL")))
dp <- Initialise(dp)
# Error: BLAS/LAPACK routine 'DPOTRFDSTEBZDSPEVXDS' gave error code -4
```

**Failing**: All MVNormal models (FULL and constrained) in R fallback mode  
**Root Cause**: `rWishart()` function receiving non-positive-definite matrix

### **Root Cause Analysis**

1. **Numerical Stability**: The prior parameters or posterior parameters create matrices that are not positive definite
2. **Parameter Computation**: The `PosteriorParameters.mvnormal()` function may compute invalid covariance matrices
3. **Prior Specification**: Default prior parameters may not be appropriate for all scenarios

### **Likely Locations**

- `R/mvnormal_normal_wishart.R:451` - `rWishart(n, post_parameters$nu_n, post_parameters$t_n)`
- `PosteriorParameters.mvnormal()` - Computation of `post_parameters$t_n`
- `MvnormalCreate()` - Prior parameter setup

### **Investigation Steps**

#### **Phase 1: Isolate Numerical Issues**
1. **Check Input Data Properties**:
   ```r
   # Analyze input data characteristics
   analyze_input_data_properties()
   check_covariance_matrix_properties()
   ```

2. **Examine Prior Parameters**:
   ```r
   # Check if prior parameters are causing issues
   examine_prior_parameters()
   test_different_prior_specifications()
   ```

3. **Debug Parameter Computation**:
   ```r
   # Step through posterior parameter computation
   debug_posterior_parameters()
   check_matrix_positive_definiteness()
   ```

#### **Phase 2: Identify Numerical Problems**
1. **Test Different Scenarios**:
   ```r
   # Test with different data characteristics
   test_various_data_sizes()
   test_different_covariance_structures()
   ```

2. **Check Matrix Conditioning**:
   ```r
   # Analyze matrix conditioning and numerical stability
   check_matrix_conditioning()
   test_regularization_approaches()
   ```

#### **Phase 3: Fix Numerical Issues**
1. **Improve Parameter Computation**:
   - Add numerical stability checks
   - Implement matrix regularization
   - Improve prior parameter defaults

2. **Alternative Sampling Methods**:
   - Consider alternative to `rWishart` for constrained models
   - Implement direct parameter sampling for constrained models

### **Implementation Plan**

#### **Step 1: Numerical Diagnosis**
- [ ] Create test script to analyze input data properties
- [ ] Examine prior parameter computation and defaults
- [ ] Debug posterior parameter computation step-by-step
- [ ] Identify specific matrices causing positive definiteness issues

#### **Step 2: Parameter Analysis**
- [ ] Test different prior specifications
- [ ] Analyze matrix conditioning throughout computation
- [ ] Identify numerical stability bottlenecks
- [ ] Document parameter ranges that cause issues

#### **Step 3: Numerical Fixes**
- [ ] Implement matrix regularization techniques
- [ ] Improve prior parameter defaults
- [ ] Add numerical stability checks
- [ ] Consider alternative sampling methods for constrained models

#### **Step 4: Validation**
- [ ] Test R fallback with various data scenarios
- [ ] Verify numerical stability improvements
- [ ] Ensure consistency between R and C++ when both work

---

## Systematic Troubleshooting Protocol

### **Priority Order**

1. **HIGH PRIORITY**: Fix C++ implementation (enables high-performance constrained models)
2. **MEDIUM PRIORITY**: Fix R fallback (ensures compatibility and debugging capability)

### **Comprehensive Testing Strategy**

#### **Test Matrix for All Scenarios**
```r
# Test matrix: [Model] x [Implementation] x [Data Size] x [Data Properties]
models <- c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
implementations <- c("cpp_enabled", "r_fallback")
data_sizes <- c(10, 50, 100, 500)
data_properties <- c("well_conditioned", "ill_conditioned", "high_dimensional")
```

#### **Regression Testing**
```r
# Ensure fixes don't break existing functionality
test_existing_functionality()
test_performance_regression()
test_result_consistency()
```

### **Success Criteria**

#### **Primary Goals**
- [ ] All constrained models work with C++ enabled
- [ ] R fallback works reliably as backup
- [ ] No performance regression in existing functionality

#### **Secondary Goals**
- [ ] Comprehensive test suite for all scenarios
- [ ] Documentation of numerical stability considerations
- [ ] Performance benchmarks for all implementations

### **Risk Assessment**

#### **High Risk Areas**
- **C++ Matrix Operations**: Complex to debug, may require C++ expertise
- **Numerical Stability**: May require deep understanding of Wishart distribution properties

#### **Mitigation Strategies**
- Create comprehensive test suite before making changes
- Document all changes for potential rollback
- Test extensively with various data scenarios
- Consider consulting numerical computation experts for Wishart issues

---

## Diagnostic Tools and Scripts

### **Immediate Next Steps**

1. **Create Diagnostic Scripts**:
   ```r
   # debug_scripts/diagnose_cpp_issues.R
   # debug_scripts/diagnose_numerical_issues.R
   # debug_scripts/comprehensive_model_testing.R
   ```

2. **Set Up Testing Infrastructure**:
   ```r
   # Create automated testing for all model/implementation combinations
   # Set up performance benchmarking
   # Create regression test suite
   ```

3. **Documentation**:
   ```r
   # Document current behavior and expected behavior
   # Create troubleshooting guide
   # Document fix progress and results
   ```

This detailed analysis provides a systematic approach to resolving the remaining issues and achieving full constrained covariance functionality in the Dirichlet Process package.