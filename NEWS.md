# dirichletprocess 0.4.2

## Major Enhancements

### Complete C++ Implementation
* **High-Performance Backend**: Comprehensive C++ implementations for all major distributions providing significant speedups
* **Automatic Fallback**: Robust R implementation fallback ensures reliability and compatibility
* **100% Manual MCMC Coverage**: Complete C++ support for 6 distributions (Normal, Exponential, Beta, Weibull, MVNormal, MVNormal2)
* **Unified Interface**: `CppMCMCRunner` provides consistent interface across all C++ implementations

### Advanced MCMC Algorithms
* **Neal's Algorithm 4**: Complete conjugate MCMC implementation in C++
* **Neal's Algorithm 8**: Full non-conjugate MCMC with Metropolis-Hastings sampling
* **Temperature Control**: Advanced features for research applications
* **Convergence Diagnostics**: Enhanced monitoring and validation tools

### Comprehensive Covariance Models
* **Full Matrix Support**: Unrestricted covariance matrices for multivariate normal distributions
* **Constrained Models**: Complete implementation of EII, VII, EEI, VEI, EVI, VVI covariance structures
* **Dimension-Aware Handling**: Robust parameter management for high-dimensional data
* **Performance Optimized**: Efficient algorithms for large covariance matrices

### Hierarchical Models
* **Complete Hierarchical Support**: 3 distributions (Hierarchical Beta, MVNormal, MVNormal2)  
* **Full MCMC C++**: High-performance hierarchical sampling algorithms
* **Global Parameter Management**: Sophisticated parameter sharing across groups
* **Stick-Breaking Implementation**: Advanced hierarchical process components

## Performance Improvements
* **Substantial Speedups**: C++ implementations provide significant performance gains for large datasets
* **Memory Efficiency**: Optimized memory usage for high-dimensional problems
* **Scalable Architecture**: Enhanced performance for complex models
* **Benchmarking Framework**: Comprehensive performance validation tools

## Technical Improvements
* **Enhanced Parameter Handling**: Robust validation and error checking
* **Better Error Messages**: Improved user feedback and diagnostics  
* **Comprehensive Testing**: Extensive R/C++ consistency validation
* **Production Ready**: Stable, well-tested implementations ready for research use

## Testing and Validation
* **Comprehensive Test Suite**: Extensive validation framework ensuring R/C++ consistency
* **Edge Case Testing**: Boundary conditions and error handling validation
* **Performance Benchmarks**: Systematic validation of C++ speedups
* **Integration Testing**: Complete workflow validation across all components

## Bug Fixes and Maintenance
* Fixed parameter handling edge cases in multivariate normal distributions
* Improved memory management in C++ implementations  
* Enhanced error handling and user feedback
* Updated documentation and examples

## Breaking Changes
* None - full backward compatibility maintained

# dirichletprocess 0.4.0.9000

* Added PriorFunction and PriorClusters to draw from the base measure. 
* Fixed a bug in the likelihood calculation (#21) by Filippo Fiocchi.
* Added hierarchical print statement method. 

# dirichletprocess 0.4.0

* Hierarchical Normal Models added by Giovanni Sighinolfi
* Added Giovanni Sighinolfi as a contributor. 
* Added params chain to Hidden Markov Models
* Updated the vignette for hierarchical normal models. 

# dirichletprocess 0.3.1

* Fixed matrix class checking for R 4.0.0
* Corrected typos in vignette
* Added a parameter for the number of initial clusters in DirichletProcessMvnormal.
* Various refactoring.


# dirichletprocess 0.3.0

* Added Hidden Markov models.
* Fixed bug in PosteriorClusters and PosteriorFunction.
* Added in new Beta mixture model for avoiding boundary. 
* Fixed a bug in ChangeObservations when using more than one dimension. 
* Added in Burn, Print and Diagnostic Plots.

# dirichletprocess 0.2.2

* Added a likelihood variable for the `dirichletprocess` class that is calculate with each fit iteration. 
* Added option to change how many Metropolis-Hasting steps are used in each iteration. 
* Added a likelihood calculation with each iteration. 
* Added and refactored some tests. 
* Updated vignette.
* Additional options to univariate plotting. 
* Added Kees Mulder as a contributor. 

# dirichletprocess 0.2.1

* Added AppVeyor, Travis-CI and codecov.io badges.
* Added penalised log-likelihood step for posterior cluster parameter inference.
* Added exponential mixture model `DirichletProcessExponential`. 
* Updated `plot`. Multivariate Gaussian models can now be plotted.
* Various bug fixes.
* Updated description.
 

# dirichletprocess 0.2.0

* First public release.
* Added a `NEWS.md` file to track changes to the package.



