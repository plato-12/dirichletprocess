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

  # Draw Wishart matrices with error handling
  sig <- tryCatch({
    simplify2array(lapply(seq_len(n),
                         function(x) {
                           wishart_draw <- rWishart(1, priorParameters$nu0, solve(priorParameters$phi0))[,,1]
                           # Handle potential numerical issues
                           if (any(is.na(wishart_draw)) || any(is.infinite(wishart_draw))) {
                             return(diag(ncol(priorParameters$phi0)))
                           }
                           tryCatch({
                             solve(wishart_draw)
                           }, error = function(e) {
                             # If singular, return identity matrix
                             diag(ncol(priorParameters$phi0))
                           })
                         }))
  }, error = function(e) {
    # Fallback to identity matrices
    array(diag(ncol(priorParameters$phi0)), dim = c(ncol(priorParameters$phi0), ncol(priorParameters$phi0), n))
  })

  # Draw multivariate normal values with error handling
  mu <- tryCatch({
    simplify2array(lapply(seq_len(n),
                          function(x) {
                            draw <- mvtnorm::rmvnorm(1, priorParameters$mu0, priorParameters$sigma0)
                            # Handle potential NA values
                            if (any(is.na(draw))) {
                              return(as.numeric(priorParameters$mu0))
                            }
                            return(draw)
                          }))
  }, error = function(e) {
    # Fallback to prior mean
    array(rep(as.numeric(priorParameters$mu0), n), dim = c(1, length(priorParameters$mu0), n))
  })

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

#' @export
#' @rdname MhParameterProposal
MhParameterProposal.mvnormal2 <- function(mdObj, old_params) {
  
  priorParameters <- mdObj$priorParameters
  new_params <- old_params
  
  # Extract current values
  old_mu <- if (is.array(old_params[[1]]) && length(dim(old_params[[1]])) == 3) {
    old_params[[1]][, , 1]
  } else {
    old_params[[1]]
  }
  
  old_sig <- if (is.array(old_params[[2]]) && length(dim(old_params[[2]])) == 3) {
    old_params[[2]][, , 1]
  } else {
    old_params[[2]]
  }
  
  # Propose new mu using multivariate normal proposal
  new_mu <- tryCatch({
    mvtnorm::rmvnorm(1, old_mu, 0.1 * old_sig)
  }, error = function(e) {
    # If covariance matrix is singular, use identity matrix
    mvtnorm::rmvnorm(1, old_mu, 0.1 * diag(length(old_mu)))
  })
  
  # Handle NA values
  if (any(is.na(new_mu))) {
    new_mu <- old_mu
  }
  
  # Propose new sig using Wishart proposal (keep current for now)
  new_sig <- old_sig
  
  # Handle potential issues with covariance matrix
  if (any(is.na(new_sig)) || any(is.infinite(new_sig))) {
    new_sig <- diag(ncol(old_sig))
  }
  
  # Return in proper format
  new_params[[1]] <- array(new_mu, dim = c(1, length(new_mu), 1))
  new_params[[2]] <- array(new_sig, dim = c(nrow(new_sig), ncol(new_sig), 1))
  
  return(new_params)
}
