# Debug Scripts for MVNormal C++ Implementation

This directory contains debugging and testing scripts used during development and troubleshooting.

## Files

### `test_minimal_debug.R`
- **Purpose**: Minimal test to isolate C++ function crashes
- **Use**: Step-by-step testing of individual C++ functions
- **When to use**: When R sessions are crashing and you need to identify which function is problematic

### `test_safe_mvnormal.R`
- **Purpose**: Safe testing version with extensive error handling
- **Use**: Comprehensive testing without losing R session due to crashes
- **When to use**: After C++ changes to verify functions work before running official tests

## Usage

```r
# To run debugging tests
source("debug_scripts/test_minimal_debug.R")
source("debug_scripts/test_safe_mvnormal.R")
```

## Note

These are temporary debugging tools. The main test suite is in `tests/testthat/test-mvnormal-cpp-comprehensive.R`.