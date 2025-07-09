# Save as analysis/class_structure.R
library(dirichletprocess)

document_s3_structure <- function() {
  sink("analysis/code_documentation/s3_structure.md")

  cat("# S3 Class Structure in dirichletprocess\n\n")

  cat("## Main S3 Classes\n\n")

  cat("### dirichletprocess\n")
  cat("Base class for all Dirichlet process objects\n\n")

  cat("#### Subclasses based on mixture distribution:\n")
  cat("- `normal`: Gaussian mixture\n")
  cat("- `beta`: Beta mixture\n")
  cat("- `mvnormal`: Multivariate normal mixture\n")
  cat("- `mvnormal2`: Alternative multivariate normal implementation\n")
  cat("- `weibull`: Weibull mixture\n")
  cat("- `exponential`: Exponential mixture\n\n")

  cat("#### Subclasses based on implementation type:\n")
  cat("- `conjugate`: For conjugate prior-likelihood pairs\n")
  cat("- `nonconjugate`: For non-conjugate prior-likelihood pairs\n")
  cat("- `hierarchical`: For hierarchical Dirichlet processes\n")
  cat("- `markov`: For hidden Markov models\n\n")

  cat("### MixingDistribution\n")
  cat("Class for mixture distribution objects\n\n")

  cat("## Method Dispatch\n\n")
  cat("The package uses S3 method dispatch extensively. Key generics include:\n\n")

  generics <- c("Fit", "ClusterComponentUpdate", "ClusterParameterUpdate",
                "PriorDraw", "PosteriorDraw", "UpdateAlpha", "Likelihood")

  for (g in generics) {
    cat(paste0("### ", g, "\n"))
    # List methods for this generic
    methods <- as.character(methods(g))
    methods <- methods[grepl(paste0("^", g, "\\."), methods)]

    if (length(methods) > 0) {
      cat("Methods:\n")
      for (m in methods) {
        cat(paste0("- `", m, "`\n"))
      }
    } else {
      cat("No specific methods found\n")
    }
    cat("\n")
  }

  sink()

  cat("S3 class structure documentation complete.\n")
}

document_s3_structure()
