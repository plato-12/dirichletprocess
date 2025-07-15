# Covariance Models Implementation Plan

## 🎯 **Objective**
Fix the missing S3 method dispatch and backend integration for covariance models (E, V, EII, VII, EEI, VEI, EVI, VVI) in the dirichletprocess package.

## 🔍 **Current Status**

### ✅ **What's Working**
- `MvnormalCreate()` function validates all covariance models
- `extractCovarianceParams()` and `reconstructCovarianceMatrix()` helper functions exist
- Model-specific parameter handling implemented
- Documentation complete

### ❌ **Critical Issue**
```
Error: no applicable method for 'Initialise' applied to an object of class 
"c('list', 'dirichletprocess', 'MixingDistribution', 'NonHierarchical')"
```

**Root Cause**: When `covModel` parameter is used, the resulting object doesn't have proper S3 method dispatch for the dirichletprocess framework.

## 📋 **Step-by-Step Implementation Plan**

### **Phase 1: Investigate Current Architecture**

#### Step 1.1: Analyze Existing S3 Methods
- **File**: `R/initialise.R`
- **Task**: Check what `Initialise` methods exist
- **Command**: `grep -n "Initialise\." R/initialise.R`
- **Expected**: Find `Initialise.mvnormal` method

#### Step 1.2: Examine Class Structure
- **File**: `R/mvnormal_normal_wishart.R`
- **Task**: Check how `MvnormalCreate()` sets object classes
- **Current**: `class(mdObj) <- c("mvnormal", "MixingDistribution", "NonHierarchical")`
- **Issue**: Missing covariance model-specific classes

#### Step 1.3: Check Method Dispatch Pattern
- **Files**: `R/normal_inverse_gamma.R`, `R/exponential_gamma.R`
- **Task**: Compare how other distributions handle S3 methods
- **Goal**: Understand the pattern for method dispatch

### **Phase 2: Fix Class Hierarchy**

#### Step 2.1: Update MvnormalCreate Class Assignment
- **File**: `R/mvnormal_normal_wishart.R`
- **Location**: Line ~110 in `MvnormalCreate()`
- **Current Code**:
```r
class(mdObj) <- c("mvnormal", "MixingDistribution", "NonHierarchical")
```
- **New Code**:
```r
# Add covariance model-specific class
if (priorParameters$covModel != "FULL") {
  class(mdObj) <- c(paste0("mvnormal.", priorParameters$covModel), 
                    "mvnormal", "MixingDistribution", "NonHierarchical")
} else {
  class(mdObj) <- c("mvnormal", "MixingDistribution", "NonHierarchical")
}
```

#### Step 2.2: Alternative Approach - Use Attributes
- **Option**: Instead of class hierarchy, use attributes
- **Code**:
```r
class(mdObj) <- c("mvnormal", "MixingDistribution", "NonHierarchical")
attr(mdObj, "covModel") <- priorParameters$covModel
```

### **Phase 3: Implement Missing S3 Methods**

#### Step 3.1: Add Initialise Methods for Constrained Models
- **File**: `R/initialise.R`
- **Task**: Add methods for each covariance model
- **Template**:
```r
#' @export
#' @rdname Initialise
Initialise.mvnormal.E <- function(mdObj, dpObj) {
  # Call base mvnormal initialise with covariance model handling
  return(Initialise.mvnormal(mdObj, dpObj))
}

#' @export
#' @rdname Initialise
Initialise.mvnormal.V <- function(mdObj, dpObj) {
  return(Initialise.mvnormal(mdObj, dpObj))
}

# ... repeat for EII, VII, EEI, VEI, EVI, VVI
```

#### Step 3.2: Update Base Initialise.mvnormal Method
- **File**: `R/initialise.R`
- **Task**: Ensure base method handles covariance models
- **Check**: Look for covariance model-specific initialization logic
- **Add**: Model-specific parameter initialization if needed

### **Phase 4: Fix Other Required S3 Methods**

#### Step 4.1: Update PosteriorDraw Methods
- **File**: `R/mvnormal_normal_wishart.R`
- **Task**: Add covariance model-specific methods
- **Pattern**:
```r
#' @export
#' @rdname PosteriorDraw
PosteriorDraw.mvnormal.EII <- function(mdObj, x, n = 1, ...) {
  # Handle EII-specific posterior draw logic
  return(PosteriorDraw.mvnormal(mdObj, x, n, ...))
}
```

#### Step 4.2: Update PriorDraw Methods
- **File**: `R/mvnormal_normal_wishart.R`
- **Task**: Add covariance model-specific methods
- **Note**: May need different parameter handling for each model

#### Step 4.3: Update Likelihood Methods
- **File**: `R/mvnormal_normal_wishart.R`
- **Task**: Ensure likelihood calculation handles all models
- **Check**: Verify `reconstructCovarianceMatrix()` is used correctly

### **Phase 5: Parameter Handling Fixes**

#### Step 5.1: Fix Parameter Dimensions
- **Issue**: Different covariance models have different parameter counts
- **Models**:
  - `E, V`: 1 parameter
  - `EII, VII`: 1 parameter
  - `EEI, EVI, VVI`: d parameters
  - `VEI`: d+1 parameters
  - `FULL`: d(d+1)/2 parameters

#### Step 5.2: Update Parameter Extraction
- **File**: `R/mvnormal_normal_wishart.R`
- **Function**: `extractCovarianceParams()`
- **Task**: Ensure all models extract correct number of parameters

#### Step 5.3: Update Parameter Reconstruction
- **File**: `R/mvnormal_normal_wishart.R`
- **Function**: `reconstructCovarianceMatrix()`
- **Task**: Verify all models reconstruct covariance correctly

### **Phase 6: C++ Integration**

#### Step 6.1: Check C++ Covariance Support
- **Files**: `src/MVNormalDistribution.cpp`, `inst/include/MVNormalDistribution.h`
- **Task**: Verify C++ implementation supports covariance models
- **Check**: Look for `covModel` parameter handling

#### Step 6.2: Update C++ Wrappers
- **File**: `R/cpp_mvnormal_wrappers.R`
- **Task**: Ensure C++ functions handle covariance models
- **Check**: Parameter conversion between R and C++

### **Phase 7: Testing and Validation**

#### Step 7.1: Create Unit Tests
- **File**: `tests/testthat/test_covariance_models.R`
- **Task**: Test each model individually
- **Tests**:
  - Object creation
  - Parameter validation
  - MCMC fitting
  - Likelihood calculation

#### Step 7.2: Integration Testing
- **File**: `debug_scripts/test_covariance_integration.R`
- **Task**: Test full pipeline for each model
- **Include**: Benchmark-style testing

#### Step 7.3: Performance Testing
- **Task**: Compare performance across models
- **Metrics**: Speed, memory usage, clustering quality

## 🔧 **Implementation Priority**

### **High Priority (Must Fix)**
1. **Fix class hierarchy** (Step 2.1)
2. **Add Initialise methods** (Step 3.1)
3. **Test basic functionality** (Step 7.1)

### **Medium Priority (Important)**
4. **Fix parameter handling** (Step 5.1-5.3)
5. **Update other S3 methods** (Step 4.1-4.3)
6. **Integration testing** (Step 7.2)

### **Low Priority (Enhancement)**
7. **C++ integration** (Step 6.1-6.2)
8. **Performance testing** (Step 7.3)

## 📝 **Implementation Template**

### **File Structure**
```
R/
├── mvnormal_normal_wishart.R     # Main implementation
├── initialise.R                  # Add Initialise methods
├── cpp_mvnormal_wrappers.R       # C++ integration
tests/testthat/
├── test_covariance_models.R      # Unit tests
debug_scripts/
├── test_covariance_integration.R # Integration tests
├── debug_covariance_models.R     # Debug utilities
```

### **Code Template for Each Model**
```r
# In R/initialise.R
#' @export
#' @rdname Initialise
Initialise.mvnormal.{MODEL} <- function(mdObj, dpObj) {
  # Model-specific initialization if needed
  return(Initialise.mvnormal(mdObj, dpObj))
}

# In R/mvnormal_normal_wishart.R
#' @export
#' @rdname PosteriorDraw
PosteriorDraw.mvnormal.{MODEL} <- function(mdObj, x, n = 1, ...) {
  # Model-specific posterior draw logic
  return(PosteriorDraw.mvnormal(mdObj, x, n, ...))
}
```

## 🚨 **Critical Files to Modify**

1. **`R/mvnormal_normal_wishart.R`** - Update class assignment in `MvnormalCreate()`
2. **`R/initialise.R`** - Add `Initialise.mvnormal.{MODEL}` methods
3. **`R/cpp_mvnormal_wrappers.R`** - Update C++ integration
4. **`tests/testthat/test_covariance_models.R`** - Create comprehensive tests

## 🎯 **Success Criteria**

- [ ] All covariance models create objects successfully
- [ ] `DirichletProcessMvnormal()` works with all models
- [ ] MCMC fitting completes without errors
- [ ] Parameter extraction/reconstruction works correctly
- [ ] C++ integration functions properly
- [ ] Unit tests pass for all models
- [ ] Benchmark script can use all models

## 🔄 **Next Steps**

1. **Start with Phase 1** - Investigate current architecture
2. **Implement Phase 2** - Fix class hierarchy (quickest fix)
3. **Test immediately** - Run covariance models test after each phase
4. **Iterate** - Fix issues as they arise
5. **Document** - Update this plan as implementation progresses

## 📊 **Estimated Timeline**

- **Phase 1-2**: 2-4 hours (investigation + class fix)
- **Phase 3-4**: 4-6 hours (S3 methods implementation)
- **Phase 5**: 2-3 hours (parameter handling)
- **Phase 6**: 3-4 hours (C++ integration)
- **Phase 7**: 2-3 hours (testing)

**Total Estimated Time**: 13-20 hours

## 🎉 **End Goal**

All covariance models working in the benchmark script:
```r
BENCHMARK_CONFIG <- list(
  covariance_models = c("FULL", "E", "V", "EII", "VII", "EEI", "VEI", "EVI", "VVI")
)
```

And successful test output:
```
=== RECOMMENDATION ===
✓ ALL COVARIANCE MODELS WORKING!
You can enable all models in the benchmark script.
```