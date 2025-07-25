# Debug hierarchical beta fitting issue
library(dirichletprocess)

# Disable C++ to isolate R implementation issue
set_use_cpp(FALSE)

# Create simple test data
group1_data <- rbeta(5, 2, 5)
group2_data <- rbeta(5, 5, 2)
hierarchical_data <- list(group1_data, group2_data)

# Create hierarchical beta DP
dp <- DirichletProcessHierarchicalBeta(hierarchical_data, maxY = 1)

cat("Starting hierarchical beta debugging...\n")
cat("Created hierarchical DP with", length(dp$indDP), "individual DPs\n")

# Try one MCMC step manually
tryCatch({
  cat("Attempting ClusterComponentUpdate on hierarchical object...\n")
  dp_updated <- ClusterComponentUpdate(dp)
  cat("ClusterComponentUpdate succeeded\n")
}, error = function(e) {
  cat("ClusterComponentUpdate failed with error:", e$message, "\n")
  cat("Error details:\n")
  print(e)
})

# Try on individual DP
tryCatch({
  cat("Attempting ClusterComponentUpdate on individual DP 1...\n")
  dp1_updated <- ClusterComponentUpdate(dp$indDP[[1]])
  cat("Individual DP ClusterComponentUpdate succeeded\n")
}, error = function(e) {
  cat("Individual DP ClusterComponentUpdate failed with error:", e$message, "\n")
  cat("Error details:\n")
  print(e)
  
  # More detailed debugging
  cat("Individual DP structure:\n")
  cat("- Class:", class(dp$indDP[[1]]), "\n")
  cat("- Has aux:", "aux" %in% names(dp$indDP[[1]]), "\n")
  if ("aux" %in% names(dp$indDP[[1]])) {
    cat("- aux[[1]] class:", class(dp$indDP[[1]]$aux[[1]]), "\n")
    cat("- aux[[1]] names:", names(dp$indDP[[1]]$aux[[1]]), "\n")
    if ("mu" %in% names(dp$indDP[[1]]$aux[[1]])) {
      cat("- aux[[1]]$mu exists\n")
    } else {
      cat("- aux[[1]]$mu MISSING\n")
    }
    if ("nu" %in% names(dp$indDP[[1]]$aux[[1]])) {
      cat("- aux[[1]]$nu exists\n")
    } else {
      cat("- aux[[1]]$nu MISSING\n")
    }
  }
})

cat("Debugging complete\n")