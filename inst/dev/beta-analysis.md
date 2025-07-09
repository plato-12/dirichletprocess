# Beta Distribution Analysis for C++ Implementation

## Key Components from R Implementation

### 1. Parameterization
- **mu**: Mean parameter (0 < mu < maxT)
- **tau**: Precision parameter (tau > 0)
- Transformation to standard Beta parameters:
  - a = (mu * tau) / maxT
  - b = (1 - mu/maxT) * tau

### 2. Prior Distribution
- mu ~ Uniform(0, maxT)
- nu ~ InverseGamma(priorParameters[1], priorParameters[2])
  - Note: nu = 1/tau (stored as nu in the code)

### 3. Likelihood Function
```r
y = (1/maxT) * dbeta(x/maxT, a, b)
```

### 4. Metropolis-Hastings Components
- **Proposal**: Random walk with step sizes for mu and nu
- **Prior Density**: Product of uniform and inverse gamma densities
- **Acceptance Ratio**: Standard MH ratio

### 5. Key Functions to Implement
1. `likelihood()` - Beta likelihood calculation
2. `priorDraw()` - Draw from prior distribution
3. `priorDensity()` - Calculate prior density
4. `mhParameterProposal()` - Generate proposal for MH
5. `posteriorDraw()` - MH sampler for posterior
6. `updatePriorParameters()` - Update hyperparameters

### 6. Special Considerations
- maxT parameter for bounded support
- Numerical stability in likelihood calculations
- Proper handling of boundary conditions in proposals