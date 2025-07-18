# Load and prepare ZIP digit recognition dataset
# Original dataset: https://web.stanford.edu/~hastie/ElemStatLearn/datasets/zip.train.gz
# This is the dataset mentioned in GitHub issue #18

load_zip_data <- function(data_path = "datasets/zip.train", 
                          max_samples = NULL, 
                          digits = NULL,
                          remove_label = TRUE) {
  
  # Read the dataset
  cat("Loading ZIP digit recognition dataset...\n")
  zip_data <- read.table(data_path, header = FALSE)
  
  # First column is the digit label (0-9)
  labels <- zip_data[, 1]
  features <- zip_data[, -1]  # 256 features (16x16 pixel intensities)
  
  cat("Original dataset dimensions:", nrow(zip_data), "x", ncol(zip_data), "\n")
  cat("Feature dimensions:", nrow(features), "x", ncol(features), "\n")
  cat("Digit distribution:\n")
  print(table(labels))
  
  # Filter by specific digits if requested
  if (!is.null(digits)) {
    keep_idx <- labels %in% digits
    labels <- labels[keep_idx]
    features <- features[keep_idx, ]
    cat("Filtered to digits:", paste(digits, collapse = ", "), "\n")
    cat("Filtered dimensions:", nrow(features), "x", ncol(features), "\n")
  }
  
  # Subsample if requested
  if (!is.null(max_samples) && max_samples < nrow(features)) {
    sample_idx <- sample(nrow(features), max_samples)
    labels <- labels[sample_idx]
    features <- features[sample_idx, ]
    cat("Subsampled to", max_samples, "samples\n")
  }
  
  # Convert to matrix for dirichletprocess
  features_matrix <- as.matrix(features)
  
  # Return data structure
  result <- list(
    data = features_matrix,
    labels = labels,
    n_samples = nrow(features_matrix),
    n_features = ncol(features_matrix),
    digit_counts = table(labels)
  )
  
  cat("Final dataset ready:", result$n_samples, "samples x", result$n_features, "features\n")
  return(result)
}

# Convenience function for high-dimensional benchmarking
prepare_benchmark_data <- function(dimensions = c(2, 10, 50, 256), 
                                   sample_sizes = c(100, 500, 1000),
                                   digits = c(0, 1, 2, 3, 4)) {
  
  # Load full dataset
  full_data <- load_zip_data(digits = digits)
  
  # Prepare data subsets for different dimensions and sample sizes
  benchmark_datasets <- list()
  
  for (dim in dimensions) {
    for (n in sample_sizes) {
      # Skip if not enough samples
      if (n > full_data$n_samples) {
        cat("Warning: Requested", n, "samples but only", full_data$n_samples, "available\n")
        next
      }
      
      # Skip if more dimensions than available
      if (dim > full_data$n_features) {
        cat("Warning: Requested", dim, "dimensions but only", full_data$n_features, "available\n")
        next
      }
      
      # Create subset
      sample_idx <- sample(full_data$n_samples, n)
      feature_idx <- if (dim == full_data$n_features) {
        1:full_data$n_features
      } else {
        sample(full_data$n_features, dim)
      }
      
      subset_data <- full_data$data[sample_idx, feature_idx, drop = FALSE]
      subset_labels <- full_data$labels[sample_idx]
      
      # Store with descriptive name
      name <- paste0("d", dim, "_n", n)
      benchmark_datasets[[name]] <- list(
        data = subset_data,
        labels = subset_labels,
        dimensions = dim,
        sample_size = n,
        name = name
      )
    }
  }
  
  cat("Prepared", length(benchmark_datasets), "benchmark datasets\n")
  return(benchmark_datasets)
}

# Quick test of the data loading
if (FALSE) {
  # Test data loading
  test_data <- load_zip_data(max_samples = 100, digits = c(0, 1, 2))
  print(str(test_data))
  
  # Test benchmark preparation
  benchmark_data <- prepare_benchmark_data(
    dimensions = c(2, 10), 
    sample_sizes = c(50, 100),
    digits = c(0, 1, 2)
  )
  print(names(benchmark_data))
}