# Test File Debugging Framework

## Problem Analysis

### Issue Summary
After reorganizing test files by moving C++ consistency tests to `tests/testthat/cpp-consistency/` subdirectory, basic functionality tests are failing due to missing dependencies.

### Root Cause
1. **Missing Original Helper**: The original `tests/testthat/helper-testing.R` was deleted and moved to `cpp-consistency/helper-testing.R`
2. **MCMC Chain Storage Issues**: Tests expect MCMC chains to be stored but they appear to be empty (length 0)
3. **Dependency Conflicts**: The cpp-consistency helper has different function implementations that may conflict

### Specific Test Failures
From `test_dirichlet_process.R`:
- `length(dpobj$clusterParametersChain) not equal to 10` (Expected: 10, Actual: 0)
- `length(dpobj$priorParametersChain) not equal to 10` (Expected: 10, Actual: 0)

## Debugging Framework

### Phase 1: Identify Missing Dependencies
1. **Check Original Helper Functions**: Restore original helper-testing.R that basic tests depend on
2. **Identify MCMC Chain Issues**: Determine why chains aren't being stored
3. **Test Individual Components**: Validate each distribution works independently

### Phase 2: Systematic Error Resolution
1. **Minimal Reproduction**: Create simple test cases to isolate issues
2. **Component-by-Component**: Fix issues without breaking existing functionality
3. **Integration Testing**: Ensure fixes work across all test files

### Phase 3: Validation and Cleanup
1. **Full Test Suite**: Run all basic tests to ensure no regressions
2. **C++ Consistency**: Verify cpp-consistency tests still work
3. **Performance Validation**: Ensure no performance degradation

## Specific Fixes Needed

### 1. Restore Original Helper Functions
- Need to recover the original helper-testing.R that basic tests depend on
- Current cpp-consistency/helper-testing.R has different functionality

### 2. MCMC Chain Storage Issues
- Tests expect `clusterParametersChain` and `priorParametersChain` to be populated
- These appear to be empty (length 0) after Fit() calls
- Need to check if this is a C++ vs R implementation issue

### 3. Test Environment Isolation
- Ensure cpp-consistency tests don't interfere with basic tests
- Separate helper functions for different test purposes

## Implementation Strategy

### Step 1: Emergency Restoration
```bash
# Restore original helper-testing.R from git
git checkout HEAD~1 -- tests/testthat/helper-testing.R
```

### Step 2: Isolate Test Dependencies
- Keep original helper-testing.R for basic tests
- Ensure cpp-consistency has its own isolated helper functions
- No name conflicts between helper functions

### Step 3: Systematic Validation
1. Run basic tests individually to identify specific issues
2. Check C++ vs R implementation differences
3. Fix chain storage issues if they exist

### Step 4: Integration Testing
- Ensure all basic tests pass
- Verify cpp-consistency tests still work
- Full test suite validation

## Test Categories and Expected Behavior

### Basic Functionality Tests (Should Pass)
- Core distribution implementations (normal, exponential, beta, etc.)
- MCMC sampling and chain storage
- Object creation and manipulation
- S3 method dispatch

### C++ Consistency Tests (Separate Directory)
- R vs C++ implementation comparison
- Performance benchmarking
- Statistical consistency validation
- Advanced C++ feature testing

## Monitoring and Validation

### Success Criteria
1. All basic functionality tests pass without errors
2. C++ consistency tests remain functional in their subdirectory
3. No test environment contamination between basic and cpp-consistency tests
4. MCMC chains are properly stored and accessible

### Risk Mitigation
1. Keep backup of working cpp-consistency tests
2. Test each fix incrementally
3. Maintain git history for rollback capability
4. Document all changes for future reference

## Next Steps
1. Restore original helper-testing.R
2. Run diagnostic tests to identify specific MCMC chain issues
3. Fix chain storage problems systematically
4. Validate all tests pass without interference