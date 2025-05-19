# dirichletprocess Core Workflow

## 1. Object Creation
The main entry point is typically through one of these constructors:
- `DirichletProcessGaussian()`
- `DirichletProcessBeta()`
- `DirichletProcessMvnormal()`
- `DirichletProcessWeibull()`

These functions all call `DirichletProcessCreate()` which sets up the basic object structure.

## 2. Initialization
The `Initialise()` function prepares the object for fitting:
- For conjugate mixtures: `Initialise.conjugate()`
- For non-conjugate mixtures: `Initialise.nonconjugate()`

## 3. Fitting
The `Fit()` function runs the MCMC sampling:
- For standard DPs: `Fit.default()`
- For hierarchical DPs: `Fit.hierarchical()`

## 4. Core MCMC Components
Each iteration of the MCMC process involves:
- `ClusterComponentUpdate()`: Updates cluster assignments
- `ClusterParameterUpdate()`: Updates cluster parameters
- `UpdateAlpha()`: Updates concentration parameter

