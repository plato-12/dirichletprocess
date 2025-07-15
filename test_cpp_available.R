library(dirichletprocess)
set_use_cpp(TRUE)

cpp_available <- function() {
  if (!using_cpp()) return(FALSE)
  
  tryCatch({
    test_md <- MvnormalCreate()
    test_data <- matrix(c(0, 0), nrow = 1)
    test_theta <- list(mu = c(0, 0), sig = diag(2))
    result <- Likelihood(test_md, test_data, test_theta)
    return(is.numeric(result) && length(result) > 0)
  }, error = function(e) {
    cat("Error in cpp_available:", e$message, "\n")
    return(FALSE)
  })
}

cat("cpp_available result:", cpp_available(), "\n")