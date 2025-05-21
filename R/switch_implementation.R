#' Set Implementation Preference
#'
#' Switch between R and C++ implementations of core functions
#'
#' @param use_cpp Logical indicating whether to use C++ implementations
#' @return Previous setting (invisibly)
#' @export
#' @examples
#' old_setting <- set_use_cpp(TRUE)
#' # Operations will now use C++ where available
#' set_use_cpp(old_setting) # Restore previous setting
set_use_cpp <- function(use_cpp) {
  old <- getOption("dirichletprocess.use_cpp", FALSE)
  options(dirichletprocess.use_cpp = use_cpp)
  invisible(old)
}

#' Check Implementation Preference
#'
#' Check whether R or C++ implementations are being used
#'
#' @return Logical indicating if C++ implementations are being used
#' @export
using_cpp <- function() {
  getOption("dirichletprocess.use_cpp", FALSE)
}

#' Get Appropriate Implementation
#'
#' @param func_name Name of the function to get implementation for
#' @return Function implementation (R or C++ based on setting)
#' @keywords internal
get_implementation <- function(func_name) {
  if (using_cpp()) {
    cpp_env <- get("cpp_implementations", envir = parent.env(environment()))
    if (exists(func_name, envir = cpp_env)) {
      return(get(func_name, envir = cpp_env))
    }
    warning("No C++ implementation available for ", func_name,
            ", falling back to R implementation")
  }
  # Return the R implementation (which is the default)
  get(func_name, envir = parent.env(environment()))
}
