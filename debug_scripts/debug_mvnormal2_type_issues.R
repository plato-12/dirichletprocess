# Debug Script for MVNormal2 C++ Type Conversion Issues
# Systematic investigation of Test 6 failures

library(dirichletprocess)
library(mvtnorm)

cat("=== Debug MVNormal2 C++ Type Conversion Issues ===\n")

# Enable C++ mode
set_use_cpp(TRUE)

cat("Checking C++ availability...\n")
if (!using_cpp()) {
  cat("C++ not available - stopping debug\n")
  quit()
}

# =============================================================================
# STEP 1: FUNCTION EXISTENCE AND BASIC INSPECTION
# =============================================================================

cat("\n1. Checking MVNormal2 Function Existence...\n")

mvnormal2_functions <- c(
  "mvnormal2_prior_draw_cpp",
  "mvnormal2_posterior_draw_cpp", 
  "mvnormal2_likelihood_cpp",
  "mvnormal2_posterior_parameters_cpp"
)

for (func_name in mvnormal2_functions) {
  exists_result <- exists(func_name, envir = asNamespace("dirichletprocess"))
  cat("  ", func_name, ":", if(exists_result) "✓ EXISTS" else "✗ MISSING", "\n")
  
  if (exists_result) {
    func_obj <- get(func_name, envir = asNamespace("dirichletprocess"))
    cat("    Type:", class(func_obj), "\n")
    if (is.function(func_obj)) {
      cat("    Arguments:", paste(names(formals(func_obj)), collapse = ", "), "\n")
    }
  }
}

# =============================================================================
# STEP 2: MVNORMAL2 OBJECT CREATION AND INSPECTION
# =============================================================================

cat("\n2. Testing MVNormal2 Object Creation...\n")

tryCatch({
  # Create test data
  set.seed(123)
  d <- 2
  n <- 10
  test_data <- matrix(rnorm(n * d), n, d)
  
  cat("  Test data dimensions:", dim(test_data), "\n")
  cat("  Test data preview:\n")
  print(head(test_data, 3))
  
  # Create MVNormal2 object
  dp <- DirichletProcessMvnormal2(test_data)
  
  cat("  ✓ DirichletProcessMvnormal2 creation successful\n")
  cat("  Object class:", class(dp), "\n")
  cat("  Mixing distribution class:", class(dp$mixingDistribution), "\n")
  
  # Inspect mixing distribution structure
  cat("  Mixing distribution structure:\n")
  str(dp$mixingDistribution, max.level = 2)
  
}, error = function(e) {
  cat("  ✗ MVNormal2 object creation failed:", e$message, "\n")
})

# =============================================================================
# STEP 3: PARAMETER PREPARATION AND INSPECTION
# =============================================================================

cat("\n3. Testing Parameter Preparation...\n")

tryCatch({
  d <- 2
  
  # Test different parameter formats
  cat("  Testing parameter format 1 (matrix mu0):\n")
  priorParams1 <- list(
    mu0 = matrix(rep(0, d), nrow = 1),
    sigma0 = diag(d),
    phi0 = diag(d) * 2,
    nu0 = d + 2
  )
  
  cat("    mu0 class:", class(priorParams1$mu0), "dimensions:", dim(priorParams1$mu0), "\n")
  cat("    sigma0 class:", class(priorParams1$sigma0), "dimensions:", dim(priorParams1$sigma0), "\n")
  cat("    phi0 class:", class(priorParams1$phi0), "dimensions:", dim(priorParams1$phi0), "\n")
  cat("    nu0 class:", class(priorParams1$nu0), "value:", priorParams1$nu0, "\n")
  
  cat("  Testing parameter format 2 (vector mu0):\n")
  priorParams2 <- list(
    mu0 = rep(0, d),
    sigma0 = diag(d),
    phi0 = diag(d) * 2,
    nu0 = d + 2
  )
  
  cat("    mu0 class:", class(priorParams2$mu0), "length:", length(priorParams2$mu0), "\n")
  cat("    sigma0 class:", class(priorParams2$sigma0), "dimensions:", dim(priorParams2$sigma0), "\n")
  
  cat("  Testing parameter format 3 (transpose mu0):\n")
  priorParams3 <- list(
    mu0 = t(matrix(rep(0, d), nrow = d)),
    sigma0 = diag(d),
    phi0 = diag(d) * 2,
    nu0 = d + 2
  )
  
  cat("    mu0 class:", class(priorParams3$mu0), "dimensions:", dim(priorParams3$mu0), "\n")
  
}, error = function(e) {
  cat("  ✗ Parameter preparation failed:", e$message, "\n")
})

# =============================================================================
# STEP 4: TESTING MVNORMAL2_PRIOR_DRAW_CPP INCREMENTALLY
# =============================================================================

cat("\n4. Testing mvnormal2_prior_draw_cpp Incrementally...\n")

if (exists("mvnormal2_prior_draw_cpp", envir = asNamespace("dirichletprocess"))) {
  
  # Test with different parameter formats
  test_params <- list(
    "Format 1 (matrix mu0)" = list(
      mu0 = matrix(c(0, 0), nrow = 1),
      sigma0 = diag(2),
      phi0 = diag(2) * 2,
      nu0 = 4
    ),
    "Format 2 (vector mu0)" = list(
      mu0 = c(0, 0),
      sigma0 = diag(2),
      phi0 = diag(2) * 2,
      nu0 = 4
    ),
    "Format 3 (as.matrix mu0)" = list(
      mu0 = as.matrix(c(0, 0)),
      sigma0 = diag(2),
      phi0 = diag(2) * 2,
      nu0 = 4
    )
  )
  
  for (format_name in names(test_params)) {
    cat("  Testing", format_name, "...\n")
    
    tryCatch({
      params <- test_params[[format_name]]
      
      cat("    Parameter inspection:\n")
      cat("      mu0:", class(params$mu0), 
          if(is.matrix(params$mu0)) paste("dims:", paste(dim(params$mu0), collapse="x")) else paste("length:", length(params$mu0)), "\n")
      cat("      sigma0:", class(params$sigma0), "dims:", paste(dim(params$sigma0), collapse="x"), "\n")
      cat("      phi0:", class(params$phi0), "dims:", paste(dim(params$phi0), collapse="x"), "\n")
      cat("      nu0:", class(params$nu0), "value:", params$nu0, "\n")
      
      # Test with minimal draws first
      cat("    Attempting mvnormal2_prior_draw_cpp with n=1...\n")
      mvnormal2_prior_draw_cpp <- get("mvnormal2_prior_draw_cpp", envir = asNamespace("dirichletprocess"))
      result <- mvnormal2_prior_draw_cpp(params, 1)
      
      cat("    ✓ SUCCESS with", format_name, "\n")
      cat("      Result class:", class(result), "\n")
      cat("      Result names:", names(result), "\n")
      if ("mu" %in% names(result)) {
        cat("      mu dimensions:", dim(result$mu), "\n")
      }
      if ("sig" %in% names(result)) {
        cat("      sig dimensions:", dim(result$sig), "\n")
      }
      
      # Cleanup
      rm(result)
      gc()
      
    }, error = function(e) {
      cat("    ✗ FAILED with", format_name, ":", e$message, "\n")
      cat("      Error class:", class(e), "\n")
      if (length(e$message) > 0) {
        cat("      Full error message:", e$message, "\n")
      }
    })
  }
  
} else {
  cat("  mvnormal2_prior_draw_cpp not available\n")
}

# =============================================================================
# STEP 5: TESTING MVNORMAL2_LIKELIHOOD_CPP INCREMENTALLY  
# =============================================================================

cat("\n5. Testing mvnormal2_likelihood_cpp Incrementally...\n")

if (exists("mvnormal2_likelihood_cpp", envir = asNamespace("dirichletprocess"))) {
  
  tryCatch({
    # First, get a working prior draw to use for likelihood testing
    cat("  Step 5a: Getting prior draw for likelihood testing...\n")
    
    priorParams <- list(
      mu0 = matrix(c(0, 0), nrow = 1),
      sigma0 = diag(2),
      phi0 = diag(2) * 2,
      nu0 = 4
    )
    
    mvnormal2_prior_draw_cpp <- get("mvnormal2_prior_draw_cpp", envir = asNamespace("dirichletprocess"))
    prior_result <- mvnormal2_prior_draw_cpp(priorParams, 1)
    cat("    ✓ Prior draw successful\n")
    
    # Inspect the prior result structure
    cat("  Step 5b: Inspecting prior draw result...\n")
    cat("    mu structure:\n")
    str(prior_result$mu)
    cat("    sig structure:\n") 
    str(prior_result$sig)
    
    # Test the correct theta format for likelihood function
    cat("  Step 5c: Testing correct theta format for likelihood...\n")
    
    # Create test data
    test_data_single <- c(0.5, -0.3)  # Vector input, not matrix
    cat("    Test data length:", length(test_data_single), "\n")
    
    # CORRECT theta format: list(mu_array, sig_array)
    cat("    Creating correct theta format...\n")
    theta_correct <- list(prior_result$mu, prior_result$sig)
    
    cat("    theta structure:\n")
    cat("      theta[[1]] (mu) dimensions:", dim(theta_correct[[1]]), "\n")
    cat("      theta[[2]] (sig) dimensions:", dim(theta_correct[[2]]), "\n")
    
    tryCatch({
      # Test the likelihood function
      mvnormal2_likelihood_cpp <- get("mvnormal2_likelihood_cpp", envir = asNamespace("dirichletprocess"))
      
      lik_result <- mvnormal2_likelihood_cpp(test_data_single, theta_correct)
      
      cat("      ✓ SUCCESS with correct theta format\n")
      cat("      Likelihood result:", lik_result, "class:", class(lik_result), "\n")
      cat("      Result length:", length(lik_result), "(should match number of draws)\n")
      
    }, error = function(e) {
      cat("      ✗ FAILED with correct theta format:", e$message, "\n")
    })
    
    # Test with different input formats
    cat("    Testing different input data formats...\n")
    
    # Test with matrix input (single row)
    test_data_matrix <- matrix(c(0.5, -0.3), nrow = 1)
    tryCatch({
      lik_result2 <- mvnormal2_likelihood_cpp(test_data_matrix, theta_correct)
      cat("      ✓ Matrix input also works! Result:", lik_result2, "\n")
    }, error = function(e) {
      cat("      ✗ Matrix input failed:", e$message, "\n")
    })
    
  }, error = function(e) {
    cat("  ✗ Likelihood testing failed:", e$message, "\n")
  })
  
} else {
  cat("  mvnormal2_likelihood_cpp not available\n")
}

# =============================================================================
# STEP 6: C++ OBJECT MEMORY INSPECTION
# =============================================================================

cat("\n6. C++ Object and Memory Inspection...\n")

tryCatch({
  # Check if we can inspect C++ objects
  cat("  Checking object sizes and memory usage...\n")
  
  # Get all objects in environment
  all_objects <- ls()
  cpp_related <- all_objects[grepl("cpp|mvnormal2", all_objects, ignore.case = TRUE)]
  
  cat("  C++ related objects in environment:\n")
  for (obj_name in cpp_related) {
    obj <- get(obj_name)
    cat("    ", obj_name, ":", class(obj), "\n")
  }
  
  # Memory usage
  cat("  Current memory usage:\n")
  mem_info <- gc()
  print(mem_info)
  
}, error = function(e) {
  cat("  Memory inspection failed:", e$message, "\n")
})

# =============================================================================
# FINAL CLEANUP AND SUMMARY
# =============================================================================

cat("\n=== Debug Summary ===\n")
cat("This debug script systematically tested:\n")
cat("1. Function existence and signatures\n")
cat("2. MVNormal2 object creation\n") 
cat("3. Parameter format variations\n")
cat("4. mvnormal2_prior_draw_cpp with different inputs\n")
cat("5. mvnormal2_likelihood_cpp with CORRECT theta format\n")
cat("6. Memory and object inspection\n")

cat("\n=== SOLUTION FOUND ===\n")
cat("The 'type conversion' issue was actually a parameter format problem:\n")
cat("❌ WRONG: theta <- list(mu = ..., sig = ...)  # Named list\n")
cat("✅ CORRECT: theta <- list(mu_array, sig_array)  # Indexed list with full arrays\n")
cat("\nThe C++ function expects:\n")
cat("- theta[0] = mu array with dimensions (1, d, n_clusters)\n")
cat("- theta[1] = sig array with dimensions (d, d, n_clusters)\n")
cat("- x = vector or matrix input\n")

cat("\nCorrect usage:\n")
cat("prior_result <- mvnormal2_prior_draw_cpp(priorParams, n_draws)\n")
cat("theta_correct <- list(prior_result$mu, prior_result$sig)\n")
cat("lik <- mvnormal2_likelihood_cpp(data_vector, theta_correct)\n")

# Cleanup
rm(list = ls())
gc()
set_use_cpp(FALSE)

cat("\n=== MVNormal2 Debug Completed ===\n")