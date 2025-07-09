# S3 Class Structure in dirichletprocess

## Main S3 Classes

### dirichletprocess
Base class for all Dirichlet process objects

#### Subclasses based on mixture distribution:
- `normal`: Gaussian mixture
- `beta`: Beta mixture
- `mvnormal`: Multivariate normal mixture
- `mvnormal2`: Alternative multivariate normal implementation
- `weibull`: Weibull mixture
- `exponential`: Exponential mixture

#### Subclasses based on implementation type:
- `conjugate`: For conjugate prior-likelihood pairs
- `nonconjugate`: For non-conjugate prior-likelihood pairs
- `hierarchical`: For hierarchical Dirichlet processes
- `markov`: For hidden Markov models

### MixingDistribution
Class for mixture distribution objects

## Method Dispatch

The package uses S3 method dispatch extensively. Key generics include:

### Fit
Methods:
- `Fit.default`
- `Fit.hierarchical`
- `Fit.markov`

### ClusterComponentUpdate
Methods:
- `ClusterComponentUpdate.conjugate`
- `ClusterComponentUpdate.hierarchical`
- `ClusterComponentUpdate.nonconjugate`

### ClusterParameterUpdate
Methods:
- `ClusterParameterUpdate.conjugate`
- `ClusterParameterUpdate.nonconjugate`

### PriorDraw
Methods:
- `PriorDraw.beta`
- `PriorDraw.beta2`
- `PriorDraw.exponential`
- `PriorDraw.hierarchical`
- `PriorDraw.mvnormal`
- `PriorDraw.mvnormal2`
- `PriorDraw.normal`
- `PriorDraw.normalFixedVariance`
- `PriorDraw.weibull`

### PosteriorDraw
Methods:
- `PosteriorDraw.exponential`
- `PosteriorDraw.mvnormal`
- `PosteriorDraw.mvnormal2`
- `PosteriorDraw.nonconjugate`
- `PosteriorDraw.normal`
- `PosteriorDraw.normalFixedVariance`
- `PosteriorDraw.weibull`

### UpdateAlpha
Methods:
- `UpdateAlpha.default`
- `UpdateAlpha.hierarchical`

### Likelihood
Methods:
- `Likelihood.beta`
- `Likelihood.beta2`
- `Likelihood.exponential`
- `Likelihood.mvnormal`
- `Likelihood.mvnormal2`
- `Likelihood.normal`
- `Likelihood.normalFixedVariance`
- `Likelihood.weibull`

