#' Create a multivariate normal mixing distribution with semi conjugate prior
#'
#' @param priorParameters The prior parameters for the Multivariate Normal.
#' @export
Mvnormal2Create <- function(priorParameters) {


  if (!is.matrix(priorParameters$mu0)){
    priorParameters$mu0 <- matrix(priorParameters$mu0, nrow=1)
  }

  mdObj <- MixingDistribution("mvnormal2", priorParameters, "nonconjugate")

  return(mdObj)
}

#' @export
#' @rdname Likelihood
Likelihood.mvnormal2 <- function(mdObj, x, theta) {
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }
  # Get dimensions and handle both 2D and 3D parameter arrays
  theta1_dim <- dim(theta[[1]])
  theta2_dim <- dim(theta[[2]])
  
  # Determine number of clusters from mu array
  num_clusters <- if (length(theta1_dim) >= 3) theta1_dim[3] else 1
  
  y <- vapply(seq_len(num_clusters),
              function(i) {
                # Extract mu for cluster i
                if (length(theta1_dim) >= 3) {
                  mu_i <- theta[[1]][, , i]
                } else {
                  mu_i <- theta[[1]][, i]
                }
                
                # Extract sigma for cluster i - handle different dimensions
                if (length(theta2_dim) >= 3) {
                  # Full covariance model - 3D array
                  sigma_i <- theta[[2]][, , i]
                } else if (length(theta2_dim) == 2) {
                  # Constrained covariance models - 2D array
                  sigma_i <- theta[[2]][, i]
                } else {
                  # Single cluster case
                  sigma_i <- theta[[2]]
                }
                
                mvtnorm::dmvnorm(x, mu_i, sigma_i)
              },
              numeric(nrow(x)))

  return(y)
}

#' @export
#' @rdname PriorDraw
PriorDraw.mvnormal2 <- function(mdObj, n = 1) {

  priorParameters <- mdObj$priorParameters

  sig <- simplify2array(lapply(seq_len(n),
                               function(x) solve(rWishart(1,
                                                          priorParameters$nu0,
                                                          solve(priorParameters$phi0))[,,1])))

  mu <- simplify2array(lapply(seq_len(n),
                              function(x) mvtnorm::rmvnorm(1,
                                                           priorParameters$mu0,
                                                           priorParameters$sigma0)))

  theta <- list(mu = mu, sig = sig)
  return(theta)
}


# PosteriorDraw.mvnormal2 <- function(mdObj, x, n = 1) {
#
#   post_parameters <- PosteriorParameters(mdObj, x)
#
#   sig <- rWishart(n, post_parameters$nu_n, post_parameters$t_n)
#   mu <- simplify2array(lapply(seq_len(n), function(x) mvtnorm::rmvnorm(1, post_parameters$mu_n,
#                                                                        solve(post_parameters$kappa_n * sig[, , x]))))
#
#   return(list(mu = mu, sig = sig/post_parameters$kappa_n^2))
# }

#' @export
#' @rdname PosteriorDraw
PosteriorDraw.mvnormal2 <- function(mdObj, x, n = 1, ...) {

  if (!is.matrix(x)) {
    x <- matrix(x, ncol = length(x))
  }

  phi0 <- mdObj$priorParameters$phi0

  mu0 <- mdObj$priorParameters$mu0
  sigma0 <- mdObj$priorParameters$sigma0

  muSamples <- array(dim = c(dim(mu0), n))
  sigSamples <- array(dim = c(dim(phi0), n))

  muSamp <- matrix(rep_len(0, ncol(mu0)), ncol=ncol(mu0))

  for (i in seq_len(n)){

    nuN <- nrow(x) +  mdObj$priorParameters$nu0
    phiN <- phi0 + Reduce("+", lapply(seq_len(nrow(x)),
                                      function(j) (x[j,] - c(muSamp)) %*% t(x[j,] - c(muSamp))))

    sigSamp <- solve(rWishart(1, nuN, solve(phiN))[,,1])

    sigN <- solve(solve(sigma0) + nrow(x) * solve(sigSamp))
    muN <- sigN %*% (nrow(x)*solve(sigSamp) %*% colMeans(x) + solve(sigma0) %*% c(mu0))

    muSamp <- mvtnorm::rmvnorm(1, muN, sigN)

    muSamples[,,i] <- muSamp
    sigSamples[,,i] <- sigSamp

  }

  return(list(mu=muSamples, sig=sigSamples))
}
