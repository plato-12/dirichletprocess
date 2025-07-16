# Debug script to trace DirichletProcessMvnormal step by step

library(dirichletprocess)
devtools::load_all()

# Disable C++ mode for testing
set_use_cpp(FALSE)

# Debug function to trace step by step
debug_step_by_step <- function() {
  cat("=== Debugging DirichletProcessMvnormal Step by Step ===\n")
  
  # Create simple 1D data
  test_data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  
  # Create E model distribution
  e_model <- MvnormalCreate(list(
    mu0 = c(0),
    kappa0 = 1,
    nu = 2,
    Lambda = diag(1),
    covModel = "E"
  ))
  
  cat("Test data dimensions:", dim(test_data), "\n")
  cat("E model created successfully.\n")
  
  # Step 1: Call DirichletProcessCreate
  cat("\n1. Testing DirichletProcessCreate...\n")
  
  tryCatch({
    dpobj <- DirichletProcessCreate(test_data, e_model, c(2, 4))
    cat("DirichletProcessCreate successful.\n")
    
    # Step 2: Call Initialise directly
    cat("\n2. Testing Initialise directly...\n")
    
    tryCatch({
      dpobj_init <- Initialise(dpobj, numInitialClusters = 1)
      cat("Initialise successful.\n")
      
    }, error = function(e) {
      cat("Initialise failed:", e$message, "\n")
      cat("Error details:\n")
      print(e)
      
      # Let's trace what happens inside Initialise
      cat("\nTracing Initialise internals...\n")
      
      # Check if we can call PosteriorDraw
      cat("  Testing PosteriorDraw...\n")
      tryCatch({
        post_draw <- PosteriorDraw(dpobj$mixingDistribution, dpobj$data, 1)
        cat("  PosteriorDraw successful.\n")
        cat("  PosteriorDraw mu dimensions:", dim(post_draw$mu), "\n")
        cat("  PosteriorDraw sig dimensions:", dim(post_draw$sig), "\n")
        
      }, error = function(e2) {
        cat("  PosteriorDraw failed:", e2$message, "\n")
      })
      
      # Check if we can call PriorDraw
      cat("  Testing PriorDraw...\n")
      tryCatch({
        prior_draw <- PriorDraw(dpobj$mixingDistribution, 1)
        cat("  PriorDraw successful.\n")
        cat("  PriorDraw mu dimensions:", dim(prior_draw$mu), "\n")
        cat("  PriorDraw sig dimensions:", dim(prior_draw$sig), "\n")
        
      }, error = function(e3) {
        cat("  PriorDraw failed:", e3$message, "\n")
      })
    })
    
  }, error = function(e) {
    cat("DirichletProcessCreate failed:", e$message, "\n")
  })
  
  # Step 3: Manually reproduce DirichletProcessMvnormal logic
  cat("\n3. Manually reproducing DirichletProcessMvnormal...\n")
  
  tryCatch({
    # Convert to matrix (already is)
    y <- test_data
    if(!is.matrix(y)){
      y <- matrix(y, ncol=length(y))
    }
    
    # Use provided priors
    g0Priors <- e_model$priorParameters
    
    # Create mixing distribution object
    mdobj <- MvnormalCreate(g0Priors)
    cat("  mdobj created successfully.\n")
    
    # Create DP object
    dpobj <- DirichletProcessCreate(y, mdobj, c(2, 4))
    cat("  dpobj created successfully.\n")
    
    # Initialize with 1 cluster
    dpobj <- Initialise(dpobj, numInitialClusters = 1)
    cat("  Initialise successful.\n")
    
  }, error = function(e) {
    cat("Manual reproduction failed:", e$message, "\n")
    print(e)
  })
  
  cat("\n=== Debug Complete ===\n")
}

# Run the debug
debug_step_by_step()