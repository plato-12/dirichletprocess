#!/usr/bin/env Rscript
# CRAN Readiness Check Script for dirichletprocess package
# This script performs comprehensive checks to ensure the package is ready for CRAN submission

# Load required packages with error handling
required_packages <- c("devtools", "rcmdcheck")
optional_packages <- c("spelling", "urlchecker")

# Load required packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Required package not installed:", pkg, "\nPlease install with: install.packages('", pkg, "')", sep = ""))
  }
  library(pkg, character.only = TRUE)
}

# Check optional packages availability
spelling_available <- requireNamespace("spelling", quietly = TRUE)
urlchecker_available <- requireNamespace("urlchecker", quietly = TRUE)

# ANSI color codes for output formatting
green <- "\033[32m"
red <- "\033[31m"
yellow <- "\033[33m"
blue <- "\033[34m"
reset <- "\033[0m"

cat_colored <- function(text, color = reset) {
  cat(paste0(color, text, reset, "\n"))
}

print_header <- function(text) {
  cat_colored(paste0("\n", "=" %&% nchar(text) %&% "="), blue)
  cat_colored(text, blue)
  cat_colored(paste0("=" %&% nchar(text) %&% "="), blue)
}

print_success <- function(text) {
  cat_colored(paste0("✓ ", text), green)
}

print_warning <- function(text) {
  cat_colored(paste0("⚠ ", text), yellow)
}

print_error <- function(text) {
  cat_colored(paste0("✗ ", text), red)
}

# Helper function for string concatenation
"%&%" <- function(x, y) paste0(x, y)

# Initialize results tracking
check_results <- list()
total_errors <- 0
total_warnings <- 0
total_notes <- 0

print_header("CRAN READINESS CHECK FOR DIRICHLETPROCESS PACKAGE")

# 1. Basic Package Structure Check
print_header("1. BASIC PACKAGE STRUCTURE")

required_files <- c("DESCRIPTION", "NAMESPACE", "R/", "man/")
missing_files <- c()

for (file in required_files) {
  if (file.exists(file)) {
    print_success(paste("Required file/directory exists:", file))
  } else {
    print_error(paste("Missing required file/directory:", file))
    missing_files <- c(missing_files, file)
  }
}

if (length(missing_files) == 0) {
  print_success("All required package files present")
} else {
  total_errors <- total_errors + length(missing_files)
}

# 2. DESCRIPTION File Validation
print_header("2. DESCRIPTION FILE VALIDATION")

desc <- try(read.dcf("DESCRIPTION"), silent = TRUE)
if (inherits(desc, "try-error")) {
  print_error("Cannot read DESCRIPTION file")
  total_errors <- total_errors + 1
} else {
  desc_df <- as.data.frame(desc, stringsAsFactors = FALSE)
  
  # Check required fields
  required_fields <- c("Package", "Title", "Version", "Description", "Authors@R", "License")
  for (field in required_fields) {
    if (field %in% colnames(desc_df) && !is.na(desc_df[[field]]) && desc_df[[field]] != "") {
      print_success(paste("Required field present:", field))
    } else {
      print_error(paste("Missing or empty required field:", field))
      total_errors <- total_errors + 1
    }
  }
  
  # Check version format
  if ("Version" %in% colnames(desc_df)) {
    version <- desc_df$Version
    if (grepl("^[0-9]+\\.[0-9]+\\.[0-9]+(-[0-9]+)?$", version)) {
      print_success(paste("Version format valid:", version))
    } else {
      print_warning(paste("Version format should follow MAJOR.MINOR.PATCH pattern:", version))
      total_warnings <- total_warnings + 1
    }
  }
  
  # Check title length
  if ("Title" %in% colnames(desc_df)) {
    title_length <- nchar(desc_df$Title)
    if (title_length <= 65) {
      print_success(paste("Title length appropriate:", title_length, "characters"))
    } else {
      print_warning(paste("Title too long:", title_length, "characters (should be ≤65)"))
      total_warnings <- total_warnings + 1
    }
  }
  
  # Check for development dependencies in Imports
  if ("Imports" %in% colnames(desc_df)) {
    imports <- desc_df$Imports
    dev_packages <- c("testthat", "devtools", "roxygen2", "knitr", "rmarkdown")
    imports_list <- trimws(unlist(strsplit(imports, ",")))
    imports_packages <- gsub("\\s*\\(.*\\)", "", imports_list)
    
    dev_in_imports <- intersect(imports_packages, dev_packages)
    if (length(dev_in_imports) > 0) {
      print_warning(paste("Development packages in Imports (should be in Suggests):", paste(dev_in_imports, collapse = ", ")))
      total_warnings <- total_warnings + 1
    } else {
      print_success("No development packages in Imports")
    }
  }
}

# 3. Documentation Check
print_header("3. DOCUMENTATION CHECK")

# Check if man pages exist
man_files <- list.files("man", pattern = "\\.Rd$")
r_files <- list.files("R", pattern = "\\.R$")

if (length(man_files) > 0) {
  print_success(paste("Documentation files found:", length(man_files), "Rd files"))
} else {
  print_error("No documentation files found in man/")
  total_errors <- total_errors + 1
}

# Check for exported functions without documentation
tryCatch({
  load_all(".", quiet = TRUE)
  namespace_exports <- getNamespaceExports("dirichletprocess")
  documented_functions <- gsub("\\.Rd$", "", man_files)
  
  undocumented <- setdiff(namespace_exports, documented_functions)
  if (length(undocumented) == 0) {
    print_success("All exported functions are documented")
  } else {
    print_warning(paste("Undocumented exported functions:", paste(undocumented, collapse = ", ")))
    total_warnings <- total_warnings + 1
  }
}, error = function(e) {
  print_warning("Could not check function documentation completeness")
  total_warnings <- total_warnings + 1
})

# 4. R CMD CHECK
print_header("4. R CMD CHECK")

cat_colored("Running R CMD check (this may take several minutes)...", blue)
check_result <- try({
  rcmdcheck::rcmdcheck(".", error_on = "never", quiet = TRUE)
}, silent = TRUE)

if (inherits(check_result, "try-error")) {
  print_error("R CMD check failed to run")
  total_errors <- total_errors + 1
} else {
  # Count issues
  errors <- length(check_result$errors)
  warnings <- length(check_result$warnings)
  notes <- length(check_result$notes)
  
  total_errors <- total_errors + errors
  total_warnings <- total_warnings + warnings
  total_notes <- total_notes + notes
  
  if (errors == 0) {
    print_success("No R CMD check errors")
  } else {
    print_error(paste("R CMD check errors:", errors))
    for (error in check_result$errors) {
      cat_colored(paste("  -", error), red)
    }
  }
  
  if (warnings == 0) {
    print_success("No R CMD check warnings")
  } else {
    print_warning(paste("R CMD check warnings:", warnings))
    for (warning in check_result$warnings) {
      cat_colored(paste("  -", warning), yellow)
    }
  }
  
  if (notes == 0) {
    print_success("No R CMD check notes")
  } else {
    print_warning(paste("R CMD check notes:", notes))
    for (note in check_result$notes) {
      cat_colored(paste("  -", substr(note, 1, 100), "..."), yellow)
    }
  }
}

# 5. Test Coverage Check
print_header("5. TEST COVERAGE CHECK")

tryCatch({
  if (dir.exists("tests")) {
    test_files <- list.files("tests", pattern = "\\.R$", recursive = TRUE)
    if (length(test_files) > 0) {
      print_success(paste("Test files found:", length(test_files)))
      
      # Try to run tests
      test_result <- try({
        devtools::test(".", quiet = TRUE)
      }, silent = TRUE)
      
      if (!inherits(test_result, "try-error")) {
        print_success("Tests run successfully")
      } else {
        print_warning("Some tests may have issues")
        total_warnings <- total_warnings + 1
      }
    } else {
      print_warning("No test files found")
      total_warnings <- total_warnings + 1
    }
  } else {
    print_warning("No tests directory found")
    total_warnings <- total_warnings + 1
  }
}, error = function(e) {
  print_warning("Could not check test coverage")
  total_warnings <- total_warnings + 1
})

# 6. Spell Check
print_header("6. SPELL CHECK")

tryCatch({
  if (spelling_available) {
    spell_errors <- spelling::spell_check_package(".")
    if (nrow(spell_errors) == 0) {
      print_success("No spelling errors found")
    } else {
      print_warning(paste("Potential spelling errors found:", nrow(spell_errors)))
      if (nrow(spell_errors) <= 10) {
        for (i in 1:nrow(spell_errors)) {
          cat_colored(paste("  -", spell_errors$word[i], "in", spell_errors$file[i]), yellow)
        }
      } else {
        cat_colored(paste("  - First 10 errors shown, total:", nrow(spell_errors)), yellow)
        for (i in 1:10) {
          cat_colored(paste("  -", spell_errors$word[i], "in", spell_errors$file[i]), yellow)
        }
      }
      total_warnings <- total_warnings + 1
    }
  } else {
    print_warning("spelling package not available for spell check (install with: install.packages('spelling'))")
  }
}, error = function(e) {
  print_warning("Could not perform spell check")
  total_warnings <- total_warnings + 1
})

# 7. URL Check
print_header("7. URL VALIDATION")

tryCatch({
  if (urlchecker_available) {
    url_results <- urlchecker::url_check(".")
    if (nrow(url_results) == 0) {
      print_success("All URLs are valid")
    } else {
      broken_urls <- url_results[url_results$Status != 200, ]
      if (nrow(broken_urls) == 0) {
        print_success("All URLs are accessible")
      } else {
        print_warning(paste("Broken URLs found:", nrow(broken_urls)))
        for (i in 1:min(5, nrow(broken_urls))) {
          cat_colored(paste("  -", broken_urls$URL[i], "(Status:", broken_urls$Status[i], ")"), yellow)
        }
        total_warnings <- total_warnings + 1
      }
    }
  } else {
    print_warning("urlchecker package not available for URL validation (install with: install.packages('urlchecker'))")
  }
}, error = function(e) {
  print_warning("Could not perform URL check")
  total_warnings <- total_warnings + 1
})

# 8. License Check
print_header("8. LICENSE CHECK")

if (file.exists("LICENSE") || file.exists("LICENCE")) {
  print_success("License file found")
} else {
  # Check if license is specified in DESCRIPTION
  if (exists("desc_df") && "License" %in% colnames(desc_df)) {
    license <- desc_df$License
    if (grepl("GPL-[23]", license) || grepl("MIT", license) || grepl("Apache", license)) {
      print_success(paste("Standard license specified:", license))
    } else {
      print_warning(paste("Non-standard license specified:", license))
      total_warnings <- total_warnings + 1
    }
  } else {
    print_error("No license information found")
    total_errors <- total_errors + 1
  }
}

# 9. News/Changelog Check
print_header("9. NEWS/CHANGELOG CHECK")

news_files <- c("NEWS.md", "NEWS", "ChangeLog", "CHANGELOG.md", "CHANGELOG")
news_found <- FALSE

for (news_file in news_files) {
  if (file.exists(news_file)) {
    print_success(paste("News file found:", news_file))
    news_found <- TRUE
    break
  }
}

if (!news_found) {
  print_warning("No NEWS or CHANGELOG file found (recommended for CRAN)")
  total_warnings <- total_warnings + 1
}

# 10. C++ Code Check (if applicable)
print_header("10. C++ CODE CHECK")

if (dir.exists("src")) {
  cpp_files <- list.files("src", pattern = "\\.(cpp|cc|c)$", recursive = TRUE)
  if (length(cpp_files) > 0) {
    print_success(paste("C++ source files found:", length(cpp_files)))
    
    # Check for Makevars
    if (file.exists("src/Makevars") || file.exists("src/Makevars.win")) {
      print_success("Makevars file(s) found")
    } else {
      print_warning("No Makevars files found (may be needed for C++ compilation)")
      total_warnings <- total_warnings + 1
    }
    
    # Check DESCRIPTION for LinkingTo
    if (exists("desc_df") && "LinkingTo" %in% colnames(desc_df)) {
      print_success("LinkingTo field found in DESCRIPTION")
    } else {
      print_warning("No LinkingTo field in DESCRIPTION (needed for C++ packages)")
      total_warnings <- total_warnings + 1
    }
  }
} else {
  print_success("No C++ source code (pure R package)")
}

# FINAL SUMMARY
print_header("CRAN READINESS SUMMARY")

if (total_errors == 0 && total_warnings == 0 && total_notes == 0) {
  cat_colored("🎉 EXCELLENT! Package appears ready for CRAN submission!", green)
} else if (total_errors == 0) {
  cat_colored("✅ GOOD! Package has no critical errors, but consider addressing warnings and notes.", yellow)
} else {
  cat_colored("❌ ISSUES FOUND! Package needs fixes before CRAN submission.", red)
}

cat_colored(paste("\nSummary:"), blue)
cat_colored(paste("Errors:  ", total_errors), if (total_errors == 0) green else red)
cat_colored(paste("Warnings:", total_warnings), if (total_warnings == 0) green else yellow)
cat_colored(paste("Notes:   ", total_notes), if (total_notes == 0) green else yellow)

# Recommendations
print_header("RECOMMENDATIONS")

if (total_errors > 0) {
  cat_colored("CRITICAL: Fix all errors before submitting to CRAN", red)
}

if (total_warnings > 0) {
  cat_colored("IMPORTANT: Address warnings to improve package quality", yellow)
}

if (total_notes > 0) {
  cat_colored("CONSIDER: Review notes and address if possible", yellow)
}

cat_colored("\nNext steps:", blue)
cat_colored("1. Fix any critical errors identified above", blue)
cat_colored("2. Consider addressing warnings and notes", blue)
cat_colored("3. Run 'devtools::check_win_devel()' for Windows testing", blue)
cat_colored("4. Run 'devtools::check_rhub()' for additional platform testing", blue)
cat_colored("5. Review CRAN policies: https://cran.r-project.org/web/packages/policies.html", blue)
cat_colored("6. Submit to CRAN: https://cran.r-project.org/submit.html", blue)

# Create summary report file
report_file <- "debug_scripts/cran_readiness_report.txt"
cat("CRAN Readiness Check Report\n", file = report_file)
cat("Generated:", as.character(Sys.time()), "\n\n", file = report_file, append = TRUE)
cat("Summary:\n", file = report_file, append = TRUE)
cat("Errors:  ", total_errors, "\n", file = report_file, append = TRUE)
cat("Warnings:", total_warnings, "\n", file = report_file, append = TRUE)
cat("Notes:   ", total_notes, "\n", file = report_file, append = TRUE)

print_success(paste("Detailed report saved to:", report_file))

# Exit with appropriate code
if (total_errors > 0) {
  quit(status = 1)
} else if (total_warnings > 0) {
  quit(status = 2)
} else {
  quit(status = 0)
}