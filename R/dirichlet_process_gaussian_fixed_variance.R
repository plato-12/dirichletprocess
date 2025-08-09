#' Create a Dirichlet Mixture of the Gaussian Distribution with fixed variance.
#'
#'
#' @param y Data.
#' @param sigma The fixed variance
#' @param g0Priors Base Distribution Priors.
#' @param alphaPriors Prior parameter distributions for the alpha concentration parameter.
#' @param cpp Logical. Use C++ implementation if TRUE, R implementation if FALSE. Default is FALSE.
#' @return Dirichlet process object
#'
#' @export
DirichletProcessGaussianFixedVariance <- function(y,
                                                  sigma,
                                                  g0Priors = c(0, 1),
                                                  alphaPriors = c(2, 4),
                                                  cpp = FALSE) {

  mdobj <- GaussianFixedVarianceMixtureCreate(g0Priors, sigma)
  dpobj <- DirichletProcessCreate(y, mdobj, alphaPriors)
  dpobj <- Initialise(dpobj)
  
  # Set cpp preference for this object
  if (cpp) {
    options(dirichletprocess.use_cpp = TRUE)
  } else {
    options(dirichletprocess.use_cpp = FALSE)
  }
  
  return(dpobj)
}
