#' Create a Multivariate Normal mixing distribution with conjugate prior
#'
#' Creates a multivariate normal mixing distribution with Normal-Wishart conjugate prior.
#' The base measure is G_0(μ, Σ | μ_0, κ_0, ν, Λ) = N(μ | μ_0, Σ/κ_0) * IW(Σ | ν, Λ)
#'
#' @param priorParameters A list containing prior parameters:
#'   \describe{
#'     \item{mu0}{Prior mean vector}
#'     \item{kappa0}{Prior precision parameter for the mean}
#'     \item{nu}{Prior degrees of freedom for the covariance}
#'     \item{Lambda}{Prior scale matrix for the covariance}
#'   }
#' @return A mixing distribution object
#' @export
MvnormalCreate <- function(priorParameters) {

  # Set default parameters if missing
  if (missing(priorParameters) || is.null(priorParameters)) {
    priorParameters <- list(
      mu0 = c(0, 0),
      kappa0 = 1,
      nu = 3,
      Lambda = diag(2)
    )
  }

  # Handle the case where a vector is passed instead of a list
  if (is.numeric(priorParameters) && !is.list(priorParameters)) {
    # Assume it's the mean vector
    d <- length(priorParameters)
    priorParameters <- list(
      mu0 = priorParameters,
      kappa0 = 1,
      nu = d + 1,
      Lambda = diag(d)
    )
  }

  # Ensure mu0 is a vector
  if (is.matrix(priorParameters$mu0)) {
    priorParameters$mu0 <- as.vector(priorParameters$mu0)
  }

  # Ensure Lambda is a matrix
  if (is.vector(priorParameters$Lambda)) {
    d <- length(priorParameters$mu0)
    if (length(priorParameters$Lambda) == 1) {
      priorParameters$Lambda <- diag(d) * priorParameters$Lambda
    } else if (length(priorParameters$Lambda) == d^2) {
      priorParameters$Lambda <- matrix(priorParameters$Lambda, nrow = d)
    }
  }

  mdObj <- MixingDistribution("mvnormal", priorParameters, "conjugate")
  return(mdObj)
}

#' @export
#' @rdname Likelihood
Likelihood.mvnormal <- function(mdObj, x, theta) {
  if (using_cpp_samplers()) {
    # Use C++ implementation if available
    if (!is.matrix(x)) {
      x <- matrix(x, nrow = 1)
    }

    # Extract parameters from theta
    mu_array <- theta[[1]]
    sig_array <- theta[[2]]

    # Handle the array structure
    if (length(dim(mu_array)) == 3) {
      n_clusters <- dim(mu_array)[3]
      d <- dim(mu_array)[2]

      result <- numeric(nrow(x))

      # For now, use first cluster (this should be generalized)
      if (n_clusters > 0) {
        mu <- mu_array[1, , 1]
        sig <- sig_array[, , 1]

        result <- mvnormal_likelihood_cpp(x, mu, sig)
      }

      return(result)
    }
  }

  # Fallback to R implementation
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }

  y <- vapply(seq_len(dim(theta[[1]])[3]),
              function(i) mvtnorm::dmvnorm(x, theta[[1]][, , i], theta[[2]][, , i]),
              numeric(nrow(x)))

  return(y)
}

#' @export
#' @rdname PriorDraw
PriorDraw.mvnormal <- function(mdObj, n = 1) {
  if (using_cpp_samplers()) {
    return(mvnormal_prior_draw_cpp(mdObj$priorParameters, n))
  }

  # Original R implementation
  priorParameters <- mdObj$priorParameters

  sig <- rWishart(n, priorParameters$nu, priorParameters$Lambda)

  mu <- simplify2array(
    lapply(seq_len(n),
           function(x)
             mvtnorm::rmvnorm(1,
                              priorParameters$mu0,
                              solve(sig[, , x] * priorParameters$kappa0))
    )
  )

  theta <- list(mu = mu, sig = sig)
  return(theta)
}

#' @export
#' @rdname PosteriorParameters
PosteriorParameters.mvnormal <- function(mdObj, x) {
  if (using_cpp_samplers()) {
    return(mvnormal_posterior_parameters_cpp(mdObj$priorParameters, as.matrix(x)))
  }

  # R implementation
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }

  priorParameters <- mdObj$priorParameters
  n <- nrow(x)
  d <- ncol(x)

  # Handle empty data
  if (n == 0) {
    return(list(
      mu_n = priorParameters$mu0,
      kappa_n = priorParameters$kappa0,
      nu_n = priorParameters$nu,
      t_n = priorParameters$Lambda
    ))
  }

  # Sample statistics
  x_bar <- colMeans(x)

  # Posterior parameters
  kappa_n <- priorParameters$kappa0 + n
  mu_n <- (priorParameters$kappa0 * priorParameters$mu0 + n * x_bar) / kappa_n
  nu_n <- priorParameters$nu + n

  # Scatter matrix
  S <- matrix(0, d, d)
  if (n > 1) {
    S <- (n - 1) * cov(x)
  }

  # Updated scale matrix
  diff <- x_bar - priorParameters$mu0
  t_n <- priorParameters$Lambda + S +
    (priorParameters$kappa0 * n / kappa_n) * outer(diff, diff)

  return(list(
    mu_n = mu_n,
    kappa_n = kappa_n,
    nu_n = nu_n,
    t_n = t_n
  ))
}

#' @export
#' @rdname PosteriorDraw
PosteriorDraw.mvnormal <- function(mdObj, x, n = 1, ...) {
  if (using_cpp_samplers()) {
    return(mvnormal_posterior_draw_cpp(mdObj$priorParameters, as.matrix(x), n))
  }

  # Original R implementation
  post_parameters <- PosteriorParameters(mdObj, x)

  sig <- rWishart(n, post_parameters$nu_n, post_parameters$t_n)
  mu <- simplify2array(
    lapply(seq_len(n),
           function(x) mvtnorm::rmvnorm(1,
                                        post_parameters$mu_n,
                                        solve(post_parameters$kappa_n * sig[, , x]))
    )
  )

  return(list(mu = mu, sig = sig))
}

#' @export
#' @rdname Predictive
Predictive.mvnormal <- function(mdObj, x) {
  if (using_cpp_samplers()) {
    return(mvnormal_predictive_cpp(mdObj$priorParameters, as.matrix(x)))
  }

  # R implementation
  if (!is.matrix(x)) {
    x <- matrix(x, nrow = 1)
  }

  priorParameters <- mdObj$priorParameters
  n <- nrow(x)
  d <- ncol(x)
  result <- numeric(n)

  for (i in 1:n) {
    x_i <- matrix(x[i, ], nrow = 1)
    post_params <- PosteriorParameters(mdObj, x_i)

    # Multivariate t-distribution predictive
    det_Lambda <- det(priorParameters$Lambda)
    det_t_n <- det(post_params$t_n)

    ratio_det <- (det_Lambda / det_t_n)^(priorParameters$nu / 2)
    ratio_kappa <- (priorParameters$kappa0 / post_params$kappa_n)^(d / 2)

    # Multivariate gamma function ratio
    gamma_ratio <- 1
    for (j in 1:d) {
      gamma_ratio <- gamma_ratio *
        gamma((post_params$nu_n + 1 - j) / 2) / gamma((priorParameters$nu + 1 - j) / 2)
    }

    pi_const <- pi^(-d/2)

    result[i] <- pi_const * ratio_kappa * ratio_det * gamma_ratio
  }

  return(result)
}
