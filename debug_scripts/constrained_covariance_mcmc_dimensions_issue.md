# Constrained Covariance Models MCMC Dimensions Issue

## Current Problem

**Primary Issue**: Constrained covariance models (EII, VII, EEI, VEI, EVI, VVI) fail during MCMC with "incorrect number of dimensions" error.

**Status**: 
- ✅ FULL model: Works completely with MCMC
- ❌ Constrained models: Fail during first MCMC iteration with dimension error

## Error Analysis

### Observed Behavior
```
Testing EII covariance model:
  ✓ Initial LikelihoodDP: length 100 
  |  0%  ✗ ERROR: incorrect number of dimensions
```

### Key Facts
1. **Likelihood works**: All models pass initial `LikelihoodDP()` calls
2. **Initialization works**: DP objects create and initialize successfully
3. **MCMC fails**: Error occurs during first MCMC iteration
4. **C++ vs R**: Error likely occurs in parameter handling during MCMC updates

## Root Cause Investigation Plan

### Phase 1: Identify Exact Failure Point

#### Step 1: Test Individual MCMC Components
Create diagnostic script to test each MCMC step separately:

```r
# Test which specific MCMC component fails
test_mcmc_components <- function(model) {
  dp <- create_constrained_dp(model)
  
  # Test cluster component update
  tryCatch({
    dp_test <- ClusterComponentUpdate(dp)
    cat("✓ ClusterComponentUpdate succeeded\n")
  }, error = function(e) {
    cat("✗ ClusterComponentUpdate failed:", e$message, "\n")
  })
  
  # Test cluster parameter update
  tryCatch({
    dp_test <- ClusterParameterUpdate(dp)  
    cat("✓ ClusterParameterUpdate succeeded\n")
  }, error = function(e) {
    cat("✗ ClusterParameterUpdate failed:", e$message, "\n")
  })
  
  # Test alpha update
  tryCatch({
    dp_test <- UpdateAlpha(dp)
    cat("✓ UpdateAlpha succeeded\n")
  }, error = function(e) {
    cat("✗ UpdateAlpha failed:", e$message, "\n")
  })
}
```

#### Step 2: Compare Parameter Structures
Analyze difference between working (FULL) and failing (constrained) models:

```r
# Compare parameter array structures
compare_parameter_structures <- function() {
  x <- matrix(rnorm(20), ncol = 2)
  
  # FULL model (working)
  dp_full <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = "FULL")))
  dp_full <- Initialise(dp_full)
  
  # EII model (failing)
  dp_eii <- DirichletProcessCreate(x, MvnormalCreate(list(covModel = "EII")))
  dp_eii <- Initialise(dp_eii)
  
  cat("=== FULL MODEL ===\n")
  cat("mu dimensions:", dim(dp_full$clusterParameters$mu), "\n")
  cat("sig dimensions:", dim(dp_full$clusterParameters$sig), "\n")
  
  cat("=== EII MODEL ===\n")
  cat("mu dimensions:", dim(dp_eii$clusterParameters$mu), "\n")
  cat("sig dimensions:", dim(dp_eii$clusterParameters$sig), "\n")
}
```

### Phase 2: Trace Parameter Flow During MCMC

#### Step 3: Debug Parameter Updates
Add detailed logging to parameter update functions:

```r
# Test parameter updates with detailed logging
debug_parameter_updates <- function(model) {
  dp <- create_constrained_dp(model)
  
  cat("Before parameter update:\n")
  print_parameter_info(dp)
  
  # Try to manually call PosteriorDraw
  tryCatch({
    new_params <- PosteriorDraw(dp$mixingDistribution, dp$data, 1)
    cat("✓ PosteriorDraw succeeded\n")
    cat("New mu dimensions:", dim(new_params$mu), "\n") 
    cat("New sig dimensions:", dim(new_params$sig), "\n")
  }, error = function(e) {
    cat("✗ PosteriorDraw failed:", e$message, "\n")
  })
}
```

#### Step 4: Check C++ Function Calls
Test if C++ functions handle constrained parameters correctly:

```r
# Test C++ function compatibility
test_cpp_compatibility <- function(model) {
  x <- matrix(rnorm(10), ncol = 2)
  prior_params <- list(
    mu0 = c(0, 0),
    kappa0 = 1,
    nu = 3,
    Lambda = diag(2),
    covModel = model
  )
  
  # Test C++ posterior draw
  tryCatch({
    result <- mvnormal_posterior_draw_cpp(prior_params, x, 1)
    cat("✓ mvnormal_posterior_draw_cpp succeeded for", model, "\n")
    cat("Result mu dim:", dim(result$mu), "\n")
    cat("Result sig dim:", dim(result$sig), "\n")
  }, error = function(e) {
    cat("✗ mvnormal_posterior_draw_cpp failed for", model, ":", e$message, "\n")
  })
}
```

### Phase 3: Fix Parameter Dimension Handling

#### Hypothesis: Array Allocation Issue
The issue is likely in how parameter arrays are allocated/reshaped for constrained models during MCMC updates.

**Expected Fix Areas**:
1. **PosteriorDraw functions**: May return wrong array dimensions for constrained models
2. **Parameter storage**: May assume FULL model array structure
3. **C++ function calls**: May expect specific dimension formats

#### Implementation Strategy:
1. **Identify the failing function** using Phase 1 diagnostics
2. **Fix parameter dimension handling** in the identified function
3. **Ensure consistent array structures** between R and C++ implementations
4. **Test all constrained models** with the fix

## Immediate Action Items

### Priority 1: Run Diagnostics
1. Create and run diagnostic scripts from Phase 1
2. Identify which MCMC component fails (ClusterComponentUpdate, ClusterParameterUpdate, or UpdateAlpha)
3. Compare parameter structures between FULL and constrained models

### Priority 2: Fix Root Cause
1. Based on diagnostics, fix the identified function's parameter handling
2. Ensure constrained model parameters are properly formatted for C++ calls
3. Test the fix with one constrained model (EII)

### Priority 3: Verify Fix
1. Test all constrained models with the fix
2. Run comprehensive MCMC tests
3. Verify performance is acceptable

## Expected Outcomes

**Success**: All 7 covariance models (FULL + 6 constrained) work correctly with MCMC and C++

**Timeline**: Should be resolvable within 1-2 debugging sessions once the exact failure point is identified

**Risk**: Low - Clear error message and working FULL model provide good debugging foundation