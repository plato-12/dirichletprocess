# CRAN Readiness Check Script Usage Guide

## Overview

The `cran_readiness_check.R` script performs comprehensive checks to ensure your R package is ready for CRAN submission. It validates package structure, documentation, runs R CMD check, and provides actionable recommendations.

## Usage

### Basic Usage
```bash
# Run from package root directory
Rscript debug_scripts/cran_readiness_check.R
```

### From R Console
```r
# Set working directory to package root
setwd("path/to/your/package")
source("debug_scripts/cran_readiness_check.R")
```

## Prerequisites

### Required Packages
- `devtools` - Essential for R package development
- `rcmdcheck` - For running R CMD check programmatically

Install with:
```r
install.packages(c("devtools", "rcmdcheck"))
```

### Optional Packages (Recommended)
- `spelling` - For spell checking documentation and code
- `urlchecker` - For validating URLs in documentation

Install with:
```r
install.packages(c("spelling", "urlchecker"))
```

## What the Script Checks

### 1. Basic Package Structure ✅
- Validates presence of required files: DESCRIPTION, NAMESPACE, R/, man/
- Ensures proper package directory structure

### 2. DESCRIPTION File Validation ✅
- Checks all required fields (Package, Title, Version, Description, Authors@R, License)
- Validates version format (MAJOR.MINOR.PATCH)
- Checks title length (should be ≤65 characters)
- Identifies development packages incorrectly placed in Imports

### 3. Documentation Check ✅
- Verifies presence of .Rd files in man/
- Identifies undocumented exported functions
- Ensures all exported functions have documentation

### 4. R CMD Check ⚠️
- Runs comprehensive R CMD check
- Reports errors, warnings, and notes
- **Critical for CRAN submission**

### 5. Test Coverage Check ✅
- Identifies test files and structure
- Attempts to run test suite
- Reports test results

### 6. Spell Check (Optional) 📝
- Checks spelling in documentation and code
- Requires `spelling` package
- Reports potential spelling errors with locations

### 7. URL Validation (Optional) 🔗
- Validates all URLs in package documentation
- Requires `urlchecker` package
- Identifies broken or inaccessible URLs

### 8. License Check ✅
- Validates license specification
- Checks for LICENSE file or DESCRIPTION license field
- Ensures compliance with CRAN requirements

### 9. News/Changelog Check 📰
- Looks for NEWS.md, CHANGELOG, or similar files
- Recommended for CRAN submissions

### 10. C++ Code Check ✅
- Detects C++ source files
- Validates Makevars files
- Checks LinkingTo field in DESCRIPTION

## Output Interpretation

### Exit Codes
- **0**: No errors or warnings - Ready for CRAN
- **1**: Critical errors found - Must fix before submission
- **2**: Warnings found - Should address before submission

### Color-Coded Output
- 🟢 **Green (✓)**: Passed checks
- 🟡 **Yellow (⚠)**: Warnings - should address
- 🔴 **Red (✗)**: Errors - must fix

### Report File
The script generates `debug_scripts/cran_readiness_report.txt` with:
- Summary of all issues found
- Timestamp of check
- Detailed breakdown of errors, warnings, and notes

## Current Status for dirichletprocess Package

Based on the latest run:

### ✅ Strengths
- Complete package structure
- Valid DESCRIPTION file
- Extensive documentation (151 .Rd files)
- Comprehensive test suite (610 tests passed)
- Proper C++ integration with Makevars
- GPL-3 license properly specified
- NEWS.md file present

### ⚠️ Issues to Address

#### Critical Error (Must Fix)
1. **R CMD Check Failed**: The script couldn't complete R CMD check due to missing Rtools
   - **Solution**: Install Rtools for R 4.4.1 from https://cran.r-project.org/bin/windows/Rtools/

#### Warning (Should Fix)
1. **Undocumented Exported Functions**: Several functions lack documentation
   - Functions like `UpdateG0.cpp`, `AlphaPriorPosteriorPlot`, etc.
   - **Solution**: Add roxygen2 documentation or remove from NAMESPACE

## Recommendations

### Immediate Actions (Before CRAN Submission)
1. **Install Rtools** to enable R CMD check
2. **Document missing functions** or remove from exports
3. **Run R CMD check manually** and address all issues
4. **Install optional packages** for comprehensive checking:
   ```r
   install.packages(c("spelling", "urlchecker"))
   ```

### Additional CRAN Preparation Steps
1. **Windows Development Testing**:
   ```r
   devtools::check_win_devel()
   ```

2. **Multi-Platform Testing**:
   ```r
   devtools::check_rhub()
   ```

3. **Review CRAN Policies**:
   - https://cran.r-project.org/web/packages/policies.html

4. **Final Manual Check**:
   ```r
   devtools::check()
   ```

## Troubleshooting

### Common Issues

#### "Required package not installed"
Install missing required packages:
```r
install.packages(c("devtools", "rcmdcheck"))
```

#### "R CMD check failed to run"
- Ensure Rtools is installed (Windows)
- Check R installation integrity
- Verify package can be loaded: `library(dirichletprocess)`

#### "Undocumented exported functions"
Either:
- Add roxygen2 documentation: `#' @export`
- Remove from NAMESPACE if not meant to be exported

#### Test Failures
- Review test output for specific failures
- Fix failing tests before CRAN submission
- Consider skipping platform-specific tests

## Best Practices

1. **Run Early and Often**: Use this script throughout development
2. **Address Issues Incrementally**: Don't wait until submission to run checks
3. **Use Version Control**: Commit fixes after addressing each issue
4. **Document Changes**: Update NEWS.md with significant changes
5. **Test on Multiple Platforms**: Use devtools helper functions

## Integration with Development Workflow

### Pre-Commit Hook
Add to `.git/hooks/pre-commit`:
```bash
#!/bin/bash
Rscript debug_scripts/cran_readiness_check.R
if [ $? -eq 1 ]; then
    echo "CRAN readiness check failed with errors. Aborting commit."
    exit 1
fi
```

### CI/CD Integration
For GitHub Actions, add to `.github/workflows/`:
```yaml
- name: CRAN Readiness Check
  run: Rscript debug_scripts/cran_readiness_check.R
```

## Support

If you encounter issues with the script:
1. Check that you're running from the package root directory
2. Ensure all required packages are installed
3. Verify R and package installations are working
4. Review the generated report file for detailed diagnostics

For package-specific CRAN submission issues, consult:
- CRAN Policies: https://cran.r-project.org/web/packages/policies.html
- R-devel mailing list archives
- RStudio Package Development guide