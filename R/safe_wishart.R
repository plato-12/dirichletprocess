# Safe rWishart wrapper with regularization
# This provides numerical stability for edge cases where rWishart fails

#' Safe rWishart wrapper with regularization
#' 
#' This function provides a robust wrapper around rWishart that handles
#' numerical edge cases where the scale matrix becomes singular or 
#' near-singular, which can cause BLAS/LAPACK errors.
#' 
#' @param n integer, number of samples
#' @param nu numeric, degrees of freedom
#' @param Lambda matrix, scale matrix (should be positive definite)
#' @param max_attempts integer, maximum number of regularization attempts
#' @param lambda_reg numeric, regularization parameter
#' @return array of Wishart samples
#' @export
safe_rWishart <- function(n, nu, Lambda, max_attempts = 3, lambda_reg = 1e-6) {
  
  # Validate inputs
  if (!is.matrix(Lambda)) {
    stop("Lambda must be a matrix")
  }
  
  if (nrow(Lambda) != ncol(Lambda)) {
    stop("Lambda must be square")
  }
  
  if (nu <= 0) {
    stop("nu must be positive")
  }
  
  if (n <= 0) {
    stop("n must be positive")
  }
  
  # Try standard rWishart first
  for (attempt in 1:max_attempts) {
    tryCatch({
      # For first attempt, use original matrix
      if (attempt == 1) {
        Lambda_reg <- Lambda
      } else {
        # Apply regularization with increasing strength
        reg_strength <- lambda_reg * (attempt - 1)
        Lambda_reg <- regularize_matrix(Lambda, method = "ridge", lambda = reg_strength)
      }
      
      # Validate that Lambda_reg is positive definite
      if (!is_positive_definite(Lambda_reg)) {
        if (attempt == max_attempts) {
          stop("Matrix is not positive definite after regularization")
        }
        next  # Try next attempt with stronger regularization
      }
      
      # Try rWishart
      result <- rWishart(n, nu, Lambda_reg)
      
      # Success - return result
      return(result)
      
    }, error = function(e) {
      # If this is the last attempt, provide informative error
      if (attempt == max_attempts) {
        stop(paste("rWishart failed after", max_attempts, "attempts. Original error:", e$message))
      }
      # Otherwise, continue to next attempt
    })
  }
  
  # Should never reach here
  stop("Unexpected error in safe_rWishart")
}

#' Regularize a matrix for numerical stability
#' 
#' @param matrix matrix to regularize
#' @param method character, regularization method ("ridge", "spectral", "nearPD")
#' @param lambda numeric, regularization parameter
#' @return regularized matrix
regularize_matrix <- function(matrix, method = "ridge", lambda = 1e-6) {
  
  switch(method,
    "ridge" = {
      # Ridge regularization - add lambda to diagonal
      matrix + diag(nrow(matrix)) * lambda
    },
    "spectral" = {
      # Spectral regularization - ensure all eigenvalues >= lambda
      eigen_decomp <- eigen(matrix)
      eigenvals <- pmax(eigen_decomp$values, lambda)
      eigen_decomp$vectors %*% diag(eigenvals) %*% t(eigen_decomp$vectors)
    },
    "nearPD" = {
      # Use Matrix::nearPD for nearest positive definite matrix
      if (requireNamespace("Matrix", quietly = TRUE)) {
        as.matrix(Matrix::nearPD(matrix, corr = FALSE)$mat)
      } else {
        # Fallback to ridge if Matrix not available
        regularize_matrix(matrix, method = "ridge", lambda = lambda)
      }
    },
    stop("Unknown regularization method: ", method)
  )
}

#' Check if matrix is positive definite
#' 
#' @param matrix matrix to check
#' @param tol numeric, tolerance for eigenvalue check
#' @return logical, TRUE if positive definite
is_positive_definite <- function(matrix, tol = 1e-8) {
  # Check if matrix is square
  if (nrow(matrix) != ncol(matrix)) {
    return(FALSE)
  }
  
  # Check if matrix is symmetric (approximately)
  if (!isSymmetric(matrix, tol = tol)) {
    return(FALSE)
  }
  
  # Check eigenvalues
  eigenvals <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  
  # All eigenvalues should be positive
  return(all(eigenvals > tol))
}

#' Test matrix conditioning
#' 
#' @param matrix matrix to analyze
#' @return list with conditioning information
analyze_matrix_conditioning <- function(matrix) {
  # Compute condition number
  cond_num <- kappa(matrix)
  
  # Get eigenvalues
  eigenvals <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  
  # Determine conditioning status
  is_well_conditioned <- cond_num < 1e12
  is_positive_definite <- all(eigenvals > 1e-8)
  
  return(list(
    condition_number = cond_num,
    eigenvalues = eigenvals,
    min_eigenvalue = min(eigenvals),
    max_eigenvalue = max(eigenvals),
    is_well_conditioned = is_well_conditioned,
    is_positive_definite = is_positive_definite,
    needs_regularization = !is_well_conditioned || !is_positive_definite
  ))
}