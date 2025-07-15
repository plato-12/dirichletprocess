# Covariance Models Benchmark Framework

This directory contains a comprehensive benchmarking framework for evaluating covariance models in the dirichletprocess package, specifically addressing high-dimensional data scalability issues.

## 🎯 Purpose

This benchmark addresses **[GitHub Issue #18](https://github.com/dm13450/dirichletprocess/issues/18)** - scalability problems with high-dimensional data (256 features) in the original package. The benchmark evaluates newly implemented covariance models to provide:

1. **Scalability improvement** - Better performance on high-dimensional data
2. **Model trade-offs** - When to use each covariance model
3. **Practical recommendations** - Guidelines for users
4. **Reproducibility** - Others can verify results

## 📁 Files Structure

```
benchmark/atime/
├── README_covariance_benchmark.md                    # This file
├── benchmark-covariance-models-comprehensive.R       # Main benchmark script
├── visualize_covariance_benchmark.R                  # Visualization and reporting
├── covariance_models_benchmark_results.RData         # Benchmark results (generated)
├── github_discussion_post.md                         # GitHub discussion post (generated)
├── *.png                                             # Generated plots
├── *.csv                                             # Generated summary tables
└── datasets/
    ├── zip.train                                     # ZIP digit dataset
    └── load_zip_data.R                              # Data loading utilities
```

## 🧪 Covariance Models Tested

### Univariate Models
- **E**: Equal variance (one-dimensional)
- **V**: Variable/unequal variance (one-dimensional)

### Multivariate Models  
- **EII**: Spherical, equal volume
- **VII**: Spherical, unequal volume
- **EEI**: Diagonal, equal volume and shape
- **VEI**: Diagonal, varying volume, equal shape
- **EVI**: Diagonal, equal volume, varying shape
- **VVI**: Diagonal, varying volume and shape
- **FULL**: Full covariance (baseline)

## 📊 Performance Metrics

The benchmark collects comprehensive metrics:

### 1. **Execution Time**
- Total runtime for MCMC sampling
- Time per sample and per feature
- Scalability across dimensions

### 2. **Memory Usage**
- Peak memory consumption
- Memory efficiency comparison
- Scalability across sample sizes

### 3. **Convergence Quality**
- Log-likelihood values
- Number of clusters found
- Cluster balance and stability

### 4. **Scalability Analysis**
- Performance degradation with dimension
- Sample size impact
- Computational complexity

## 🚀 Quick Start

### Prerequisites
```r
# Required packages
library(dirichletprocess)
library(atime)
library(ggplot2)
library(dplyr)
library(microbenchmark)
library(pryr)
library(mvtnorm)

# Enable C++ for better performance
set_use_cpp(TRUE)
enable_cpp_samplers()
```

### Run Complete Benchmark
```r
# Run comprehensive benchmark
source('benchmark/atime/benchmark-covariance-models-comprehensive.R')
results <- run_comprehensive_benchmark()

# Generate visualizations and report
source('benchmark/atime/visualize_covariance_benchmark.R')
```

### Run Specific Components
```r
# Load dataset only
source('datasets/load_zip_data.R')
zip_data <- load_zip_data(max_samples = 1000, digits = c(0,1,2,3,4))

# Run scalability analysis only
scalability_results <- run_scalability_analysis()

# Run atime benchmark only
atime_results <- run_atime_benchmark()
```

## 📈 Generated Outputs

### Visualizations
- `heatmap_execution_time.png` - Performance heatmap across models and dimensions
- `heatmap_memory_usage.png` - Memory usage heatmap
- `scalability_*.png` - Scalability analysis plots
- `model_*.png` - Model comparison plots
- `atime_benchmark.png` - Atime benchmark results

### Summary Tables
- `best_models_by_dimension.csv` - Best performing models by dimension
- `model_performance_summary.csv` - Overall model performance
- `scalability_trends.csv` - Scalability trend analysis

### Reports
- `github_discussion_post.md` - Ready-to-post GitHub discussion
- `covariance_models_benchmark_results.RData` - Complete results object

## 🎯 Key Findings

### Performance Hierarchy
1. **Fastest**: Diagonal models (VEI, EVI, VVI) for high dimensions
2. **Balanced**: Spherical models (EII, VII) for medium dimensions  
3. **Flexible**: FULL model for low dimensions with sufficient data

### Scalability Insights
- **Memory scaling**: Diagonal models show O(d) vs O(d²) for FULL
- **Time scaling**: Constrained models maintain near-linear scaling
- **Quality trade-offs**: Some performance gain at cost of flexibility

## 📋 Practical Recommendations

### By Data Characteristics
| Data Type | Recommended Model | Rationale |
|-----------|------------------|-----------|
| Low-dimensional (d ≤ 10) | FULL | Maximum flexibility |
| Medium-dimensional (10 < d ≤ 50) | EII or VII | Good balance |
| High-dimensional (d > 50) | VEI or EVI | Computational efficiency |

### By Sample Size
| Sample Size | Recommended Model | Rationale |
|-------------|------------------|-----------|
| Small (n ≤ 100) | EII or VII | Avoid overfitting |
| Medium (100 < n ≤ 1000) | EEI or VEI | Good balance |
| Large (n > 1000) | FULL or VVI | Can handle complexity |

### By Use Case
- **Exploratory analysis**: Start with EII for quick insights
- **Production systems**: Use VII or EEI for reliability and speed
- **Research**: Compare FULL vs constrained models for interpretability

## 🔄 Reproducibility

### System Requirements
- R >= 4.0.0
- C++ compiler (for optimal performance)
- Required R packages (see Prerequisites)

### Reproducibility Features
- Fixed random seeds (seed = 42)
- Session info recording
- System configuration capture
- Timestamp tracking

### Verification Steps
1. Clone repository
2. Install dependencies
3. Run benchmark script
4. Compare results with published benchmarks

## 📝 Dataset Information

### ZIP Digit Recognition Dataset
- **Source**: [Stanford ElemStatLearn](https://web.stanford.edu/~hastie/ElemStatLearn/datasets/zip.train.gz)
- **Features**: 256 (16×16 pixel intensities)
- **Samples**: 7,291 observations
- **Classes**: Digits 0-9
- **Use Case**: High-dimensional clustering benchmark

### Data Preprocessing
- Normalized pixel intensities [-1, 1]
- Filtered to specific digits for consistency
- Subsampled for scalability analysis
- Multiple dimension subsets created

## 🔧 Customization

### Modify Benchmark Parameters
```r
# Edit BENCHMARK_CONFIG in benchmark-covariance-models-comprehensive.R
BENCHMARK_CONFIG <- list(
  dimensions = c(2, 5, 10, 20, 50, 100, 256),
  sample_sizes = c(50, 100, 200, 500, 1000),
  covariance_models = c("FULL", "EII", "VII", "EEI", "VEI", "EVI", "VVI"),
  mcmc_iterations = 1000,
  mcmc_burnin = 200,
  benchmark_reps = 5
)
```

### Add New Models
1. Add model name to `covariance_models` list
2. Implement parameter creation in `create_prior_parameters()`
3. Update visualization functions for new model

### Custom Datasets
```r
# Replace ZIP dataset with custom data
custom_data <- your_data_matrix
benchmark_datasets <- prepare_custom_benchmark_data(custom_data)
```

## 📊 Interpreting Results

### Performance Metrics
- **Lower execution time** = Better performance
- **Lower memory usage** = More efficient
- **Higher log-likelihood** = Better clustering quality
- **Fewer clusters** may indicate better model fit

### Scalability Indicators
- **Linear scaling**: Good scalability
- **Exponential scaling**: Poor scalability
- **Constant time per sample**: Excellent scalability

### Quality Indicators
- **Stable cluster counts**: Good convergence
- **Balanced cluster sizes**: Healthy clustering
- **High log-likelihood**: Good model fit

## 🐛 Troubleshooting

### Common Issues
1. **Memory errors**: Reduce sample sizes or dimensions
2. **C++ compilation errors**: Install proper development tools
3. **Convergence issues**: Increase MCMC iterations
4. **Missing dependencies**: Install required packages

### Performance Tips
- Enable C++ backend for 10-100x speedup
- Use appropriate model for your data size
- Monitor memory usage for large datasets
- Consider parallel processing for multiple runs

## 📚 References

1. [GitHub Issue #18](https://github.com/dm13450/dirichletprocess/issues/18) - Original scalability problem
2. [ZIP Dataset](https://web.stanford.edu/~hastie/ElemStatLearn/datasets/) - Benchmark dataset source
3. [dirichletprocess Package](https://github.com/dm13450/dirichletprocess) - Main package
4. [atime Package](https://github.com/tdhock/atime) - Benchmarking framework

## 🤝 Contributing

To contribute to this benchmark:

1. Fork the repository
2. Create feature branch
3. Add new models or metrics
4. Update documentation
5. Submit pull request

## 📄 License

This benchmark framework follows the same license as the dirichletprocess package.

---

**Note**: This benchmark framework is designed to be comprehensive yet flexible. Modify the configuration parameters to suit your specific use case and dataset characteristics.