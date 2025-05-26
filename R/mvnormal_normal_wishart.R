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
#' @rdname PosteriorDraw
PosteriorDraw.mvnormal <- function(mdObj, x, n = 1, ...) {
  if (using_cpp_samplers()) {
    return(mvnormal_posterior_draw_cpp(mdObj$priorParameters, as.matrix(x), n))
  }

  # Original R implementation continues...
  post_parameters <- PosteriorParameters(mdObj, x)

  sig <- rWishart(n, post_parameters$nu_n, post_parameters$t_n)
  mu <- simplify2array(
    lapply(seq_len(n),
           function(x) mvtnorm::rmvnorm(1,
                                        post_parameters$mu_n,
                                        solve(post_parameters$kappa_n * sig[, , x]))
    )
  )

  return(list(mu = mu, sig = sig/post_parameters$kappa_n^2))
}
