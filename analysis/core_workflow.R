# Save as analysis/core_workflow.R
library(dirichletprocess)

# Document the main object creation and fitting workflow
document_workflow <- function() {
  sink("analysis/core_workflow.md")

  cat("# dirichletprocess Core Workflow\n\n")

  cat("## 1. Object Creation\n")
  cat("The main entry point is typically through one of these constructors:\n")
  cat("- `DirichletProcessGaussian()`\n")
  cat("- `DirichletProcessBeta()`\n")
  cat("- `DirichletProcessMvnormal()`\n")
  cat("- `DirichletProcessWeibull()`\n\n")

  cat("These functions all call `DirichletProcessCreate()` which sets up the basic object structure.\n\n")

  cat("## 2. Initialization\n")
  cat("The `Initialise()` function prepares the object for fitting:\n")
  cat("- For conjugate mixtures: `Initialise.conjugate()`\n")
  cat("- For non-conjugate mixtures: `Initialise.nonconjugate()`\n\n")

  cat("## 3. Fitting\n")
  cat("The `Fit()` function runs the MCMC sampling:\n")
  cat("- For standard DPs: `Fit.default()`\n")
  cat("- For hierarchical DPs: `Fit.hierarchical()`\n\n")

  cat("## 4. Core MCMC Components\n")
  cat("Each iteration of the MCMC process involves:\n")
  cat("- `ClusterComponentUpdate()`: Updates cluster assignments\n")
  cat("- `ClusterParameterUpdate()`: Updates cluster parameters\n")
  cat("- `UpdateAlpha()`: Updates concentration parameter\n\n")

  sink()

  cat("Core workflow documentation complete.\n")
}

document_workflow()
