# R/zzz.R (package startup)

.onLoad <- function(libname, pkgname) {
  # Set default options
  options(
    dirichletprocess.use_cpp = TRUE,  # Default to C++ if available
    dirichletprocess.cpp_debug = FALSE  # For debugging
  )
}

#' Set whether to use C++ implementation
#' @export
set_use_cpp <- function(use_cpp = TRUE) {
  options(dirichletprocess.use_cpp = use_cpp)
  invisible(use_cpp)
}

#' Check if using C++ implementation
#' @export
using_cpp <- function() {
  getOption("dirichletprocess.use_cpp", FALSE)
}
