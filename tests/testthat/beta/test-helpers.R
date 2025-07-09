# Helper functions for tests

#' Generate mixture of Beta distributions
#' @param n Total sample size
#' @param params List of (a, b) pairs for each component
#' @param weights Mixture weights (default: equal)
generate_beta_mixture <- function(n, params, weights = NULL) {
  k <- length(params)
  if (is.null(weights)) weights <- rep(1/k, k)

  component_sizes <- as.integer(n * weights)
  component_sizes[k] <- n - sum(component_sizes[-k])  # Adjust for rounding

  data <- numeric(n)
  idx <- 1

  for (i in 1:k) {
    size <- component_sizes[i]
    if (size > 0) {
      data[idx:(idx + size - 1)] <- rbeta(size, params[[i]][1], params[[i]][2])
      idx <- idx + size
    }
  }

  sample(data)  # Shuffle
}

#' Check MCMC convergence using Geweke diagnostic
#' @param chain MCMC chain
#' @param frac1 Fraction of chain for first window
#' @param frac2 Fraction of chain for second window
geweke_diag <- function(chain, frac1 = 0.1, frac2 = 0.5) {
  n <- length(chain)
  n1 <- floor(n * frac1)
  n2 <- floor(n * frac2)

  chain1 <- chain[1:n1]
  chain2 <- chain[(n - n2 + 1):n]

  mean1 <- mean(chain1)
  mean2 <- mean(chain2)
  var1 <- var(chain1)
  var2 <- var(chain2)

  z <- (mean1 - mean2) / sqrt(var1/n1 + var2/n2)

  list(z = z, p.value = 2 * pnorm(-abs(z)))
}

#' Calculate effective sample size
#' @param chain MCMC chain
eff_sample_size <- function(chain) {
  if (length(chain) < 10) return(NA)

  acf_vals <- acf(chain, lag.max = length(chain)/4, plot = FALSE)$acf

  # Find first negative autocorrelation
  first_neg <- which(acf_vals < 0)[1]
  if (is.na(first_neg)) first_neg <- length(acf_vals)

  # ESS = n / (1 + 2 * sum of positive autocorrelations)
  sum_acf <- sum(acf_vals[2:first_neg])
  ess <- length(chain) / (1 + 2 * sum_acf)

  return(ess)
}
