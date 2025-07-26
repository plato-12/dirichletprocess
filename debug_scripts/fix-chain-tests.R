# Fix Chain Expectation Tests
# This script identifies and fixes tests that expect clusterParametersChain or priorParametersChain
# to be populated but fail with C++ implementation

library(testthat)

# Find all test files that might have chain expectation issues
test_files <- list.files("tests/testthat", pattern = "test_.*\\.R$", full.names = TRUE)

cat("Searching for chain expectation patterns in", length(test_files), "test files...\n")

chain_patterns <- c(
  "expect_equal\\(length\\(.*clusterParametersChain\\)",
  "expect_length\\(.*clusterParametersChain\\)",
  "expect_equal\\(length\\(.*priorParametersChain\\)",
  "expect_length\\(.*priorParametersChain\\)"
)

# Function to check if file contains problematic patterns
check_file_for_chains <- function(file_path) {
  content <- readLines(file_path)
  matches <- list()
  
  for (pattern in chain_patterns) {
    matching_lines <- grep(pattern, content, value = TRUE)
    if (length(matching_lines) > 0) {
      matches[[pattern]] <- list(
        lines = matching_lines,
        line_numbers = grep(pattern, content)
      )
    }
  }
  
  if (length(matches) > 0) {
    return(list(file = file_path, matches = matches))
  }
  return(NULL)
}

# Check all test files
problematic_files <- list()
for (file in test_files) {
  result <- check_file_for_chains(file)
  if (!is.null(result)) {
    problematic_files[[file]] <- result
  }
}

# Report findings
cat("\n=== CHAIN EXPECTATION ANALYSIS ===\n")
if (length(problematic_files) == 0) {
  cat("No problematic chain expectations found.\n")
} else {
  cat("Found", length(problematic_files), "files with chain expectation issues:\n\n")
  
  for (file_info in problematic_files) {
    cat("FILE:", basename(file_info$file), "\n")
    for (pattern in names(file_info$matches)) {
      cat("  Pattern:", pattern, "\n")
      for (i in seq_along(file_info$matches[[pattern]]$lines)) {
        line_num <- file_info$matches[[pattern]]$line_numbers[i]
        line_content <- file_info$matches[[pattern]]$lines[i]
        cat("    Line", line_num, ":", trimws(line_content), "\n")
      }
    }
    cat("\n")
  }
}

cat("=== ANALYSIS COMPLETE ===\n")