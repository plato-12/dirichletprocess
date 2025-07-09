# Save as analysis/create_master_doc.R
create_master_doc <- function() {
  sink("analysis/dirichletprocess_analysis.md")

  cat("# dirichletprocess Package Analysis\n\n")

  # Add structure overview
  cat("## Package Structure\n\n")
  cat(readLines("analysis/structure_overview.txt"), sep = "\n")
  cat("\n\n")

  # Add core workflow
  cat("## Core Workflow\n\n")
  cat(readLines("analysis/core_workflow.md"), sep = "\n")
  cat("\n\n")

  # Add S3 structure
  cat("## S3 Class Structure\n\n")
  cat(readLines("analysis/code_documentation/s3_structure.md"), sep = "\n")
  cat("\n\n")

  # Add mathematical basis
  cat("## Mathematical Foundation\n\n")
  cat(readLines("analysis/code_documentation/mathematical_basis.md"), sep = "\n")
  cat("\n\n")

  # Add bottleneck analysis
  cat("## Performance Analysis\n\n")
  cat(readLines("analysis/profiling_results/bottleneck_analysis.md"), sep = "\n")
  cat("\n\n")

  # Add optimization roadmap
  cat("## Optimization Plan\n\n")
  cat(readLines("analysis/optimization_roadmap.md"), sep = "\n")

  sink()

  cat("Master documentation created at analysis/dirichletprocess_analysis.md\n")
}

create_master_doc()
