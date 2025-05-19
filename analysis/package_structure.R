# analysis/package_structure.R
library(fs)
library(dplyr)
library(purrr)

# --- Main Package Analysis ---
cat("Current working directory:", getwd(), "\n")

# Get all R files
# This assumes there is a directory named "R" in your current working directory.
r_files <- fs::dir_ls("R", regexp = "\\.R$", recurse = TRUE, fail = FALSE) # Added fail = FALSE

cat("Found", length(r_files), "R files in R/ directory.\n")

# Ensure the 'analysis' directory exists before writing files
# This will be relative to getwd()
analysis_dir <- "analysis" # Definition for the main analysis part
if (!fs::dir_exists(analysis_dir)) {
  fs::dir_create(analysis_dir)
  cat("Created '", analysis_dir, "' directory in ", getwd(), ".\n")
}

# Create a dataframe with file information
if (length(r_files) > 0) {
  file_info <- tibble(
    file_path = r_files,
    file_name = fs::path_file(r_files)
  ) %>%
    mutate(
      # Fetch file information (like size and modification time) using fs::file_info
      file_details = fs::file_info(file_path),
      size = file_details$size,                   # Extract size from file_details
      last_modified = file_details$modification_time, # Use modification_time from fs::file_info
      file_details = NULL # Remove the temporary helper column
    )

  # Extract function names from each file
  file_info <- file_info %>%
    mutate(functions = map(file_path, function(path) {
      if (!fs::file_exists(path)) return(character(0))

      lines <- readLines(path, warn = FALSE)
      func_lines <- lines[grepl("^[A-Za-z0-9_\\.]+ *<- *function\\(", lines)]
      gsub(" *<- *function.*$", "", func_lines)
    }))

  # Prepare data for CSV: first mutate to create functions_count, then select
  file_listing_data <- file_info %>%
    mutate(functions_count = purrr::map_int(functions, length)) %>%
    select(file_name, size, last_modified, functions_count)

  # Save the full structure information
  rds_path <- fs::path(analysis_dir, "file_structure.rds")
  saveRDS(file_info, rds_path)
  cat("File structure RDS saved to:", fs::path_real(rds_path), "\n")

  # Write the prepared data to CSV
  csv_path <- fs::path(analysis_dir, "file_listing.csv")
  write.csv(file_listing_data, csv_path, row.names = FALSE)
  cat("File listing CSV saved to:", fs::path_real(csv_path), "\n")

  # Create a simple report
  report_path <- fs::path(analysis_dir, "structure_overview.txt")
  cat("Package Structure Overview\n",
      "=======================\n",
      "Total R files: ", nrow(file_info), "\n",
      "Total functions: ", sum(purrr::map_int(file_info$functions, length)), "\n",
      file = report_path)
  cat("Structure overview report saved to:", fs::path_real(report_path), "\n")

  print("Package structure analysis complete. Results saved in 'analysis' directory.")

} else {
  print("No R files found in the 'R' directory. Main analysis skipped.")

  report_path <- fs::path(analysis_dir, "structure_overview.txt") # analysis_dir should be defined from top
  cat("Package Structure Overview\n",
      "=======================\n",
      "Total R files: 0\n",
      "Total functions: 0\n",
      file = report_path)
  cat("Empty structure overview report saved to:", fs::path_real(report_path), "\n")

  empty_listing_df <- tibble(
    file_name = character(),
    size = numeric(),
    last_modified = as.POSIXct(character()),
    functions_count = integer()
  )
  csv_path <- fs::path(analysis_dir, "file_listing.csv") # analysis_dir should be defined from top
  write.csv(empty_listing_df, csv_path, row.names = FALSE)
  cat("Empty file listing CSV saved to:", fs::path_real(csv_path), "\n")

  print("Empty report and file listing generated in 'analysis' directory.")
}

# --- Addendum: Create Directory Tree ---
cat("\n--- Directory Tree Generation ---\n")

# Define analysis_dir here to ensure it's available for this specific block,
# even if the block is run in isolation or if the environment might have changed.
# This mirrors the definition from earlier in the script.
analysis_dir <- "analysis"

if (requireNamespace("data.tree", quietly = TRUE)) {
  library(data.tree) # Ensure it's loaded for this block

  cat("Package 'data.tree' is available. Attempting to generate directory tree...\n")
  cat("Current working directory for tree generation:", getwd(), "\n")

  # Ensure 'analysis' directory exists (it should from earlier script parts, but double check)
  # This uses the 'analysis_dir' defined just above for this block.
  if (!fs::dir_exists(analysis_dir)) {
    fs::dir_create(analysis_dir)
    cat("Created '", analysis_dir, "' directory in ", getwd(), " (should have been created earlier if main script ran).\n")
  }

  # Get all files from the current working directory, then filter.
  all_files_for_tree <- fs::dir_ls(path = ".", recurse = TRUE, type = "file", fail = FALSE)

  # Filter for specific file types relevant to a package structure or common project files
  package_files_pattern <- "(\\.R|\\.Rmd|\\.cpp|\\.h|\\.hpp|\\.Rd|DESCRIPTION|NAMESPACE|LICENSE|NEWS\\.md|README\\.md|\\.txt|\\.json|\\.yaml|\\.yml)$"
  package_files <- grep(package_files_pattern, all_files_for_tree, value = TRUE, ignore.case = TRUE)

  cat("Found", length(package_files), "files matching pattern for directory tree.\n")

  if (length(package_files) > 0) {
    # Create a data.frame with pathString and other attributes like size
    file_attribs_for_tree <- tibble(
      pathString = fs::path_norm(fs::path_rel(package_files, start = getwd())),
      size_bytes = fs::file_size(package_files)
    ) %>% mutate(
      size_readable = format(size_bytes, standard = "SI", units = "auto")
    )

    package_tree <- tryCatch({
      data.tree::as.Node(file_attribs_for_tree, pathDelimiter = "/")
    }, error = function(e) {
      cat("Error creating tree structure with data.tree: ", conditionMessage(e), "\n")
      return(NULL)
    })

    if (!is.null(package_tree) && package_tree$totalCount > 0) {
      package_tree$Do(function(node) {
        node$type <- ifelse(node$isLeaf, "File", "Directory")
        if (is.null(node$size_readable)) node$size_readable <- ""
      }, traversal = "post-order")

      output_file_path_tree <- fs::path(analysis_dir, "directory_tree.txt")

      sink(output_file_path_tree)
      cat("Package Directory Structure (Root: ", getwd(), ")\n")
      cat("Showing files matching pattern: ", package_files_pattern, "\n\n")
      print(package_tree, "type", "size_readable", limit = 300)
      sink()

      cat("Directory tree report saved to:", output_file_path_tree, "\n")
      if(fs::file_exists(output_file_path_tree)){
        cat("SUCCESS: Directory tree file found at:", fs::path_real(output_file_path_tree), "\n")
        cat("File size:", format(fs::file_size(output_file_path_tree), standard = "SI", units = "auto"), "\n")
      } else {
        cat("ERROR: Directory tree file NOT found after attempting to write to:", output_file_path_tree, "\n")
      }

    } else {
      cat("No files found to create a directory tree, or the resulting tree was empty.\n")
      output_file_path_tree <- fs::path(analysis_dir, "directory_tree.txt")
      fs::file_write("Directory tree not generated: No matching files or empty tree.",
                     path = output_file_path_tree)
      cat("Placeholder message saved to:", fs::path_real(output_file_path_tree), "\n")
    }
  } else {
    cat("No files matching the pattern for directory tree generation were found in", getwd(), "\n")
    output_file_path_tree <- fs::path(analysis_dir, "directory_tree.txt")
    fs::file_write("Directory tree not generated: No files matched the specified pattern.",
                   path = output_file_path_tree)
    cat("Placeholder message saved to:", fs::path_real(output_file_path_tree), "\n")
  }
} else {
  cat("Package 'data.tree' is NOT installed. Skipping directory tree generation.\n")
  cat("To generate the directory tree, please install it using: install.packages('data.tree')\n")
  output_file_path_tree <- fs::path(analysis_dir, "directory_tree.txt") # analysis_dir defined at start of this block
  if (fs::dir_exists(analysis_dir)) {
    fs::file_write("Directory tree not generated: 'data.tree' package not found/installed.",
                   path = output_file_path_tree)
    cat("Placeholder message about missing 'data.tree' saved to:", fs::path_real(output_file_path_tree), "\n")
  }
}
cat("\n--- End of script ---\n")
