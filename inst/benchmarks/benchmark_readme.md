# Exponential Distribution Benchmarking Guide

## Overview

This benchmarking suite provides comprehensive performance and memory profiling tools for comparing the R and C++ implementations of the exponential distribution in the dirichletprocess package.

## Components

### 1. **benchmark_exponential.R**
Main benchmarking script with the following functions:

- `benchmark_exponential_comprehensive()`: Tests performance across various dataset sizes, iteration counts, and cluster configurations
- `profile_exponential_memory()`: Detailed memory usage tracking during fitting
- `benchmark_exponential_components()`: Component-level performance analysis
- `run_and_report()`: Convenience function to run all benchmarks and generate a summary

### 2. **exponential_viz.R**
Visualization functions:

- `create_benchmark_dashboard()`: Creates a comprehensive visualization dashboard
- `create_performance_report()`: Generates a summary plot with key statistics
- `save_benchmark_plots()`: Saves all plots to disk

### 3. **test_exponential_performance.R**
Quick testing and validation script:

- `quick_exponential_benchmark()`: Fast benchmark for quick testing
- `profile_exponential_components()`: Detailed component profiling with memory tracking
- `test_statistical_equivalence()`: Validates that R and C++ produce equivalent results

## Usage

### Quick Benchmark

```r
library(dirichletprocess)
source("test_exponential_performance.R")

# Run quick benchmark
results <- quick_exponential_benchmark()
```

### Comprehensive Benchmark

```r
source("benchmark_exponential.R")

# Run full benchmark suite
results <- run_and_report()

# Access individual results
bench_data <- results$benchmarks
memory_data <- results$memory
component_data <- results$components
```

### Generate Visualizations

```r
source("benchmark_exponential.R")
source("exponential_viz.R")

# Run benchmarks
results <- run_and_report()

# Create visualizations
plots <- create_benchmark_dashboard(
  results$benchmarks,
  results$memory,
  results$components
)

# Save plots
save_benchmark_plots(plots, "benchmark_results/exponential")
```

### Memory Profiling

The C++ implementation includes built-in memory tracking:

```r
# Enable C++ implementation
set_use_cpp(TRUE)

# Clear previous tracking
clear_memory_tracking()

# Run some operations
dp <- DirichletProcessExponential(data)
dp <- Fit(dp, 100)

# Get memory report
mem_report <- get_memory_tracking()
print(mem_report)
```

## Interpreting Results

### Performance Metrics

1. **Speedup Factor**: Ratio of R execution time to C++ execution time
   - Values > 1 indicate C++ is faster
   - Typical range: 5-20x depending on problem size

2. **Memory Usage**: Peak memory consumption during fitting
   - C++ typically uses less memory due to more efficient data structures
   - Memory tracking shows internal C++ allocations

3. **Scaling Behavior**: How performance changes with dataset size
   - Both implementations should show similar O(n) scaling
   - C++ advantage typically increases with problem size

### Component Analysis

The benchmarks break down performance into three key components:

1. **Likelihood Calculation**: Evaluating the exponential density
   - Most frequent operation
   - Largest speedup potential

2. **Cluster Assignment**: Chinese Restaurant Process updates
   - Complex algorithmic component
   - Benefits from C++ memory efficiency

3. **Parameter Updates**: Posterior sampling for rate parameters
   - Conjugate updates are fast in both implementations
   - C++ reduces overhead

## Customization

### Adjust Benchmark Parameters

```r
# Custom benchmark configuration
results <- benchmark_exponential_comprehensive(
  n_obs_vec = c(1000, 5000, 10000),    # Dataset sizes
  n_iter_vec = c(500, 1000),           # Iteration counts
  n_clusters_vec = c(3, 5, 10),        # True cluster numbers
  n_reps = 10                          # Repetitions per scenario
)
```

### Add New Metrics

Extend the benchmarking functions to track additional metrics:

```r
# Example: Track convergence speed
track_convergence <- function(dp_obj, true_rates) {
  found_rates <- sort(as.numeric(dp_obj$clusterParameters[[1]]))
  return(mean(abs(found_rates - sort(true_rates))))
}
```

## Troubleshooting

### C++ Not Available

If C++ implementations are not available:

```r
# Check status
get_cpp_status()

# Ensure C++ is compiled
# Rebuild package with: R CMD INSTALL --preclean .
```

### Memory Tracking Issues

If memory tracking returns empty results:

1. Ensure you're using the C++ implementation: `set_use_cpp(TRUE)`
2. Clear tracking before starting: `clear_memory_tracking()`
3. Some operations may not allocate trackable memory

### Performance Anomalies

If results seem inconsistent:

1. Increase `n_reps` for more stable averages
2. Ensure no other processes are consuming CPU
3. Use `gc()` between runs to ensure clean memory state
4. Check that data generation is consistent (use same seed)

## Expected Results

Based on typical hardware (modern CPU, 8GB+ RAM):

| Dataset Size | Expected Speedup | Memory Reduction |
|--------------|------------------|------------------|
| 100          | 5-10x            | 20-30%           |
| 1,000        | 10-15x           | 30-40%           |
| 10,000       | 15-25x           | 40-50%           |

Actual results will vary based on:
- CPU architecture and speed
- Available memory and cache sizes
- Compiler optimizations
- Number of clusters in data
