#' @useDynLib dirichletprocess, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @title A flexible package for fitting Bayesian non-parametric models.
#' @name dirichletprocess
#' @description Create, fit and take posterior samples from a Dirichlet process.
#'
#'
#' @importFrom stats dbeta dbinom dgamma dnorm dt dunif dweibull dexp rWishart rbeta rgamma rnorm runif var quantile optim
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#'
#'
#' @keywords internal
"_PACKAGE"
