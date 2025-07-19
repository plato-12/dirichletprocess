---
marp: true
theme: default
paginate: true
style: |
  section {
    font-size: 32px;
    font-family: 'Arial', sans-serif;
  }
---

# Scaling Bayesian Nonparametrics: C++ Acceleration of dirichletprocess

**Speaker:** Priyanshu Tiwari  
**Institution:** IIT Kanpur  
**Event:** useR! 2025 Virtual Lightning Talk  
**Duration:** 3-5 minutes  

---

## Slide 1: Title and Overview

### **Scaling Bayesian Nonparametrics: C++ Acceleration of dirichletprocess**

**Priyanshu Tiwari**  
Indian Institute of Technology Kanpur  

**Lightning Talk - useR! 2025**

**Project Focus:** Transforming R package performance through strategic C++ implementation

---

## Slide 2: What is Dirichlet Process?

### **Dirichlet Process: The Infinite Mixture Model**

**Key Concept:** A flexible Bayesian method for clustering without pre-specifying cluster numbers

**Mathematical Foundation:**
- **Chinese Restaurant Process**: Elegant metaphor for infinite clustering
- **Nonparametric**: Model complexity grows with data
- **Bayesian**: Uncertainty quantification built-in

**Applications:**
- **Density estimation** for unknown distributions
- **Clustering** with automatic cluster discovery  
- **Hierarchical modeling** across groups

**Why Important:** Handles real-world complexity where cluster numbers are unknown

---

## Slide 3: Original R Package Capabilities

### **dirichletprocess Package: Comprehensive Bayesian Toolkit**

**Core Functionality:**
- **Multiple Distributions**: Normal, Beta, Exponential, Weibull, Multivariate Normal
- **Algorithmic Flexibility**: Conjugate and non-conjugate implementations
- **Hierarchical Models**: Multi-level clustering structures
- **Object-Oriented Design**: S3 classes for extensibility

**Unique Features:**
- **Neal (2000) Algorithms**: Gold standard MCMC implementations
- **Automatic Fallback**: Seamless R/C++ integration
- **9 Covariance Models**: Complete multivariate normal coverage
- **Educational Value**: Clear, readable implementations

**Research Impact:** Enables reproducible Bayesian nonparametric analysis

---

## Slide 4: Performance Bottlenecks

### **The Challenge: Computational Scaling**

**Critical Performance Issues:**
- **MCMC Algorithms**: Inherently sequential, difficult to vectorize
- **Loop-Heavy Operations**: Chinese Restaurant Process requires extensive iteration
- **Memory Inefficiency**: R's overhead grows exponentially with dataset size
- **Scalability Limits**: Large datasets become computationally prohibitive

**Real-World Impact:**
- **Research limitations**: Studies constrained by computation time
- **Educational barriers**: Students cannot experiment with realistic datasets
- **Adoption challenges**: Practitioners choose faster alternatives

**JSS Review Feedback:** "Performance limitations prevent wider adoption"

---

## Slide 5: C++ Implementation Architecture

### **Strategic C++ Integration: Performance Without Compromise**

**Architecture Principles:**
- **Preserve R Interface**: Zero disruption to existing users
- **Automatic Fallback**: Seamless integration with error handling
- **Modular Design**: Each algorithm component separately optimized
- **Memory Efficiency**: Pre-allocated structures, minimal copying

**Technical Stack:**
- **Rcpp + RcppArmadillo**: High-performance linear algebra
- **Neal (2000) Algorithms**: Identical mathematical implementations
- **Unified Interface**: CppMCMCRunner for all distributions
- **Advanced Features**: Temperature control, cluster operations

**Coverage Achievement:** **100% manual MCMC C++ coverage** across all major distributions

---

## Slide 6: Performance Results - The Numbers

### **Transformative Performance Improvements**

**Execution Speed Gains:**
- **Beta Distribution**: **108x faster** (up to 221x speedup)
- **Multivariate Normal**: **287x faster** (up to 711x speedup)  
- **Normal Distribution**: **29x faster** (up to 39x speedup)
- **Exponential**: **17x faster** (up to 19x speedup)
- **Weibull**: **11x faster** (up to 23x speedup)

**Memory Efficiency:**
- **Normal**: **272x less memory** usage
- **Exponential**: **278x less memory** usage
- **Multivariate Normal**: **201x less memory** usage

**Scalability Impact:**
- **Large datasets**: Analysis of 10-100x larger datasets now feasible
- **Research applications**: Minutes instead of hours for complex models

---

## Slide 7: Real-World Impact Example

### **Case Study: Multivariate Normal Clustering**

**Dataset:** 500 observations, 3-dimensional clustering

**R Implementation:**
- **Time**: 848 seconds (14.1 minutes)
- **Memory**: 0.5 GB
- **Practical limit**: ~1,000 observations

**C++ Implementation:**
- **Time**: 1.19 seconds
- **Memory**: 3 MB
- **Practical limit**: ~100,000+ observations

**Result:** **711x speedup**, **201x memory reduction**

**Research Impact:** Enables real-time interactive analysis and large-scale studies

---

## Slide 8: Technical Achievements

### **Engineering Excellence: Production-Ready Implementation**

**Algorithm Fidelity:**
- **Identical Results**: Statistical equivalence with R implementations
- **Numerical Stability**: Log-space calculations prevent overflow
- **Comprehensive Testing**: 200+ unit tests ensure correctness

**Software Quality:**
- **Automatic Detection**: Runtime C++ availability checking
- **Graceful Degradation**: Seamless fallback to R implementations
- **Cross-Platform**: Windows, macOS, Linux compatibility
- **CRAN Ready**: Meets all package standards

**Developer Experience:**
- **Preserved Interface**: No learning curve for existing users
- **Advanced Features**: Temperature control, manual stepping
- **Extensible Design**: Framework for future distributions

---

## Slide 9: Future Roadmap

### **Expanding Bayesian Nonparametric Capabilities**

**Immediate Goals:**
- **Complete Minor Distributions**: Beta2, Normal Fixed Variance
- **Hierarchical Optimization**: Performance gains for multi-level models
- **Parallelization**: OpenMP integration for multi-core systems

**Advanced Features:**
- **GPU Acceleration**: CUDA implementations for massive datasets
- **Streaming Algorithms**: Online learning for continuous data
- **Variational Inference**: Approximate methods for ultra-large datasets

**Community Impact:**
- **Educational Resources**: Workshops and tutorials
- **Research Collaborations**: Applied statistics partnerships
- **Open Source**: Continued development and maintenance

**Vision:** Making Bayesian nonparametrics accessible at any scale

---

## Slide 10: Call to Action

### **Join the Bayesian Nonparametric Revolution**

**Try dirichletprocess Today:**
- **CRAN Installation**: `install.packages("dirichletprocess")`
- **GitHub Development**: github.com/plato-12/dirichletprocess
- **Documentation**: Comprehensive vignettes and examples

**Get Involved:**
- **Feedback**: Share your use cases and performance experiences
- **Contributing**: Help expand C++ coverage to more distributions
- **Research**: Collaborate on methodological improvements

**Contact Information:**
- **Email**: priyanshut21@iitk.ac.in
- **GitHub**: @plato-12
- **Website**: plato-12.github.io

**Thank You!** Questions welcome - let's scale Bayesian analysis together!

---

## Alt-Text Descriptions for Visual Elements

### **Slide 2 Visual Description:**
*Diagram showing the Chinese Restaurant Process metaphor: customers (data points) choosing tables (clusters) with probability proportional to table occupancy, illustrating the infinite mixture concept.*

### **Slide 6 Visual Description:**
*Bar chart comparing R vs C++ execution times across five distributions, with logarithmic scale showing dramatic speedup improvements ranging from 11x to 287x average speedup.*

### **Slide 7 Visual Description:**
*Side-by-side comparison showing timing results: R implementation clock showing 14.1 minutes vs C++ implementation showing 1.19 seconds, with memory usage icons showing 500MB vs 3MB.*

### **Slide 8 Visual Description:**
*Architecture diagram showing the seamless integration between R interface and C++ backend, with automatic fallback mechanisms and unified CppMCMCRunner interface.*

---

## Presentation Notes for Accessibility

### **Speaking Guidelines:**
- **Pace**: Speak slowly and clearly for 3-5 minute duration
- **Verbalization**: Read all slide text aloud, including numbers and technical terms
- **Descriptions**: Audibly describe all visual elements and charts
- **Emphasis**: Use voice inflection to highlight key performance improvements

### **Technical Requirements:**
- **Font**: Large sans-serif fonts (32pt minimum)
- **Contrast**: High contrast text and background
- **Colors**: Colorblind-safe palette with additional visual indicators
- **Code**: Maximum 50 characters per line when showing code examples

### **Downloadable Format:**
This markdown file serves as the accessible downloadable resource, enabling screen readers and braille displays to access all presentation content including detailed performance data and technical specifications.

---

## References and Additional Resources

- **Neal, R.M. (2000)**: "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
- **Package Documentation**: Comprehensive vignettes and API reference
- **Benchmark Results**: Detailed performance analysis across all distributions
- **Source Code**: Full C++ implementation with extensive comments
- **Educational Materials**: Tutorials and example applications

---

*This presentation follows useR! 2025 accessibility guidelines including large fonts, high contrast design, colorblind-safe visuals, and comprehensive alt-text descriptions for all graphical elements.*