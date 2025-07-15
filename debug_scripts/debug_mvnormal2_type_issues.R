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
  exists_result <- exists(func_name)
  cat("  ", func_name, ":", if(exists_result) "✓ EXISTS" else "✗ MISSING", "\n")
  
  if (exists_result) {
    func_obj <- get(func_name)
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

if (exists("mvnormal2_prior_draw_cpp")) {
  
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

if (exists("mvnormal2_likelihood_cpp")) {
  
  tryCatch({
    # First, get a working prior draw to use for likelihood testing
    cat("  Step 5a: Getting prior draw for likelihood testing...\n")
    
    priorParams <- list(
      mu0 = matrix(c(0, 0), nrow = 1),
      sigma0 = diag(2),
      phi0 = diag(2) * 2,
      nu0 = 4
    )
    
    prior_result <- mvnormal2_prior_draw_cpp(priorParams, 1)
    cat("    ✓ Prior draw successful\n")
    
    # Inspect the prior result structure
    cat("  Step 5b: Inspecting prior draw result...\n")
    cat("    mu structure:\n")
    str(prior_result$mu)
    cat("    sig structure:\n") 
    str(prior_result$sig)
    
    # Test different ways to extract parameters for likelihood
    cat("  Step 5c: Testing parameter extraction methods...\n")
    
    extraction_methods <- list(
      "Method 1" = list(
        mu = prior_result$mu[1, , 1],
        sig = prior_result$sig[, , 1]
      ),
      "Method 2" = list(
        mu = as.vector(prior_result$mu[1, , 1]),
        sig = as.matrix(prior_result$sig[, , 1])
      ),
      "Method 3" = list(
        mu = c(prior_result$mu[1, , 1]),
        sig = prior_result$sig[, , 1]
      )
    )
    
    # Create test data
    test_data_single <- matrix(c(0.5, -0.3), nrow = 1)
    cat("    Test data dimensions:", dim(test_data_single), "\n")
    
    for (method_name in names(extraction_methods)) {
      cat("    Testing", method_name, "...\n")
      
      tryCatch({
        theta <- extraction_methods[[method_name]]
        
        cat("      theta$mu class:", class(theta$mu), 
            if(is.matrix(theta$mu)) paste("dims:", paste(dim(theta$mu), collapse="x")) else paste("length:", length(theta$mu)), "\n")
        cat("      theta$sig class:", class(theta$sig), "dims:", paste(dim(theta$sig), collapse="x"), "\n")
        cat("      theta$mu values:", paste(theta$mu, collapse=", "), "\n")
        
        # Test the likelihood function
        lik_result <- mvnormal2_likelihood_cpp(test_data_single, theta)
        
        cat("      ✓ SUCCESS with", method_name, "\n")
        cat("      Likelihood result:", lik_result, "class:", class(lik_result), "\n")
        
      }, error = function(e) {
        cat("      ✗ FAILED with", method_name, ":", e$message, "\n")
      })
    }
    
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
cat("5. mvnormal2_likelihood_cpp with different parameter extractions\n")
cat("6. Memory and object inspection\n")

cat("\nReview the output above to identify:\n")
cat("- Which parameter formats work vs fail\n")
cat("- Exact error messages and their contexts\n")
cat("- Where the type conversion issues occur\n")
cat("- Memory or object lifecycle problems\n")

# Cleanup
rm(list = ls())
gc()
set_use_cpp(FALSE)

cat("\n=== MVNormal2 Debug Completed ===\n")