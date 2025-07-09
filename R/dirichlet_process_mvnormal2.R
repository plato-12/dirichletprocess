#' Create a Dirichlet mixture of multivariate normal distributions with semi-conjugate prior.
#'
#'
#'
#' @param y Data
#' @param g0Priors Prior parameters for the base distribution.
#' @param alphaPriors Alpha prior parameters. See \code{\link{UpdateAlpha}}.
#' @export
DirichletProcessMvnormal2 <- function(y,
                                      g0Priors,
                                      alphaPriors = c(2, 4)) {

  if (!is.matrix(y)){
    y <- matrix(y, ncol=length(y))
  }

  if(missing(g0Priors)){
    # Fix: Ensure nu0 is large enough for the Wishart distribution
    d <- ncol(y)
    g0Priors <- list(nu0 = d + 2,  # Changed from 2 to d + 2
                     phi0 = diag(d),
                     mu0 = numeric(d),
                     sigma0 = diag(d))
  }

  # Validate nu0
  if(g0Priors$nu0 <= ncol(y) - 1) {
    stop(sprintf("nu0 must be greater than %d (dimension - 1) for valid Wishart distribution",
                 ncol(y) - 1))
  }

  mdobj <- Mvnormal2Create(g0Priors)
  dpobj <- DirichletProcessCreate(y, mdobj, alphaPriors)
  dpobj <- Initialise(dpobj)

  return(dpobj)
}
